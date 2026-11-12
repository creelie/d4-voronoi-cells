# run_all.jl -- runs the Julia re-implementations of the computer checks and prints a
# PASS/FAIL table with the numbers compared against the logs of the original runs.
#
#   julia --project=. run_all.jl [--tier quick|standard|full] [--llm24 DIR | --fetch-llm24 DIR]
#                                [--only REGEX] [--list]
#
#   quick     exact checks and the small branch and bounds (a few minutes)
#   standard  quick + every check that takes at most about ten minutes on one core (default)
#   full      standard + the long branch and bounds (hours; see README.md for each)
#
# --llm24 DIR        folder proofs/4_24 of the unpacked data set of de Laat, Leijenhorst and
#                    de Muinck Keizer (4TU.ResearchData, DOI 10.4121/74ce1c25-6fca-4680-8a36-e9c18e7e9594)
# --fetch-llm24 DIR  download the archive into DIR, check its SHA-256, unpack it, and use it
# Without either, the two checks that read that data set are reported as SKIPPED.
#
# --tight-log LOG replaces the log read by profile_levels_tight (see README.md).
#
# Each check writes its own log to logs/<name>.log; the table is also written to
# logs/SUMMARY.txt (logs/SUMMARY_only.txt with --only).  The exit code is 0 only if no check failed.
include(joinpath(@__DIR__, "src", "D4Checks.jl"))
using .D4Checks
using .D4Checks: GNbhd, ProfileLevels, Constants, Common, AncCheckExact, AncAnnular, AncExterior, CapCertificate, Radial, ThreePoint, LLM24, Combo, Covering, Cardinality
using Printf

const ROOT = @__DIR__
const DATA = joinpath(ROOT, "data")
const LOGS = joinpath(ROOT, "logs")
const TIERS = Dict("quick" => 1, "standard" => 2, "full" => 3)

function parse_args(args)
    o = Dict{String,Any}("tier" => "standard", "llm24" => nothing, "fetch" => nothing, "only" => nothing, "list" => false, "tight" => nothing)
    i = 1
    while i <= length(args)
        a = args[i]
        if a == "--tier"; o["tier"] = args[i+1]; i += 2
        elseif a == "--llm24"; o["llm24"] = args[i+1]; i += 2
        elseif a == "--fetch-llm24"; o["fetch"] = args[i+1]; i += 2
        elseif a == "--tight-log"; o["tight"] = args[i+1]; i += 2
        elseif a == "--only"; o["only"] = Regex(args[i+1]); i += 2
        elseif a == "--list"; o["list"] = true; i += 1
        else error("unknown argument $a") end
    end
    haskey(TIERS, o["tier"]) || error("--tier must be quick, standard or full")
    o
end

struct Entry
    name::String
    tier::Int
    supports::String
    f::Function          # io -> Report, or nothing when skipped
end

cert(p...) = joinpath(DATA, p...)
rlog(n) = joinpath(DATA, "logs", n)

"closed + dropped boxes of a radial log (used to sort the radial checks into tiers)"
function radial_boxes(logpath)
    n = 0
    for line in eachline(logpath)
        m = match(r"(\d+) boxes closed, (\d+) dropped", line)
        m !== nothing && (n += parse(Int, m.captures[1]) + parse(Int, m.captures[2]))
    end
    n
end

function entries(llm, tight=nothing)
    E = Entry[]
    add(name, tier, supports, f) = push!(E, Entry(name, tier, supports, f))
    add("anc_check_exact", 1, "explicit constants of the annular/exterior argument (anc/check_exact.py)", io -> AncCheckExact.run(io))
    add("anc_annular_constants", 1, "annular constants (anc/annular_constants.py)", io -> AncAnnular.run_annular_constants(io))
    add("anc_cap_comparison", 1, "cap comparison (anc/cap_comparison.py)", io -> AncAnnular.run_cap_comparison(io))
    add("anc_cap_profiles_26", 1, "the 17845 profiles at 26 centres (anc/cap_profiles_26.py)", io -> AncAnnular.run_cap_profiles_26(io))
    add("anc_exterior_test", 1, "exterior computation (anc/exterior_test.py)", io -> AncExterior.run(io))
    add("cap_certificate", 1, "thm:cap-inequality, rem:cap-boxes (cap_certificate/)", io -> CapCertificate.run(io))
    add("covering_bound", 1, "thm:covering-bound, cor:m22, tab:covering-bound", io -> Covering.run(io))
    add("twentyfour_close", 1, "lem:twentyfour-close (constants in 200-bit ball arithmetic)", io -> Constants.run_twentyfour_close(io))
    add("pushout", 1, "thm:G-pushout (constants in ball arithmetic and exact rationals)", io -> Constants.run_pushout(io))
    add("G_nbhd", 1, "thm:G-nbhd, app:G-nbhd, tab:G-nbhd-bins (lem:T-near-roots, thm:one-large-push at e0 = 1/48; 200-bit balls; first case + 29 bins)", io -> GNbhd.run(io))
    add("hexagon_rigidity", 1, "thm:hexagon-rigidity, cor:hexagon-close (exact rationals)", io -> GNbhd.run_hexagon(io))
    add("labelled_certificate_II", 1, "thm:m23-labelled, regions II_s and II_f (region I is in the full tier)",
        io -> ThreePoint.run_labelled(cert("certificates", "certificate_d8.npz"), io; skip_region_1=true))
    add("certificate_check", 2, "thm:certificate (m = 23 three-point certificate)", io -> ThreePoint.run_certificate(cert("certificates", "certificate_d8.npz"), io))
    add("labelled_certificate_full", 3, "thm:m23-labelled, all three regions (region I: about 420 000 boxes)",
        io -> ThreePoint.run_labelled(cert("certificates", "certificate_d8.npz"), io))
    add("radial_count_31", 2, "thm:count31", io -> Radial.run_count(cert("radial_certificates", "radial_31.json"), io;
        expect=Dict("pair" => (191740, 180310, 54248), "intervals" => 356, "total" => "3.0853743 (Python log)")))
    # radial_case_check.py: the single splits, the tree of thm:count30 and the levels
    for f in sort(readdir(joinpath(DATA, "radial_certificates")))
        m = match(r"^(case|density)_(\d+.*)\.json$", f)
        (m === nothing || endswith(f, "_spec.json") || occursin("spec", f)) && continue
        logname = (m.captures[1] == "case" ? "radial_case_check_" : "density_check_") * m.captures[2] * ".log"
        lp = rlog(logname)
        isfile(lp) || continue
        nb = radial_boxes(lp)
        tier = nb <= 20_000 ? 1 : nb <= 600_000 ? 2 : 3
        sup = m.captures[1] == "density" ? "prop:levels, two-point level at M = $(m.captures[2])" :
              f == "case_30.json" ? "thm:count30 (the tree; about 185 million boxes)" :
              "prop:C-radial / prop:G-radial, single split $(m.captures[2])"
        add(replace(f, ".json" => ""), tier, sup * " [$nb boxes]",
            io -> Radial.run_case(cert("radial_certificates", f), io; expect=Radial.expect_from_log(lp)))
    end
    cc(n) = cert("cardinality_certificates", n)
    add("cardinality_d8_s005", 2, "thm:kissing-stable, degree 8, t = 0.505 (115 532 boxes)",
        io -> Cardinality.run(cc("cert_d8_t0.50500.npz"), io; e1="1e-6", e2="1e-5", wmin=1e-5, expect=("24.899999963", 197, 115532)))
    add("cardinality_d8_s0065", 2, "thm:kissing-stable, degree 8, t = 0.5065 (150 225 boxes)",
        io -> Cardinality.run(cc("cert_d8_t0.50650.npz"), io; e1="1e-6", e2="1e-5", wmin=1e-5, expect=("24.974999950", 221, 150225)))
    add("cardinality_d10_s008", 3, "thm:kissing-stable, degree 10, t = 0.508 (2 642 578 boxes)",
        io -> Cardinality.run(cc("cert_d10_t0.50800.npz"), io; e1="1e-6", e2="7e-5", wmin=1e-6, expect=("24.959999949", 1479, 2642578)))
    add("cardinality_26_d10", 3, "thm:twenty-six, degree 10, t = 0.51468, N = 26 (1 287 272 boxes)",
        io -> Cardinality.run(cc("cert_d10_t0.51468.npz"), io; e1="1e-6", e2="2.9e-4", wmin=1e-5, N=26, expect=("25.799999928", 790, 1287272)))
    add("combo24_level", 2, "prop:levels at M = 24 (combined two-/three-point certificate)",
        io -> Combo.run(cert("combo", "combo24_level_d6.npz"), 6, io; case=cert("combo", "case24_all.json"), border=true, level="3.3291", expect_bound=3.253872))
    add("combo30_precheck", 2, "thm:count30, fifth case: exact positivity, thresholds and bound only (no branch and bound; not a proof)",
        io -> Combo.run(cert("combo", "combo30_d8.npz"), 8, io; mg2=2e-5, mg3=5e-6, mgm=2e-6, bb=false, expect_bound=3.087203))
    add("combo30", 3, "thm:count30, fifth case (8.3 million triple boxes)",
        io -> Combo.run(cert("combo", "combo30_d8.npz"), 8, io; mg2=2e-5, mg3=5e-6, mgm=2e-6, expect_bound=3.087203))
    add("combo26_level", 3, "prop:levels at M = 26 (7.5 million triple boxes)",
        io -> Combo.run(cert("combo", "combo26_level_d8.npz"), 8, io; case=cert("combo", "cases", "case26_all.json"), level="3.3352", expect_bound=3.329050))
    add("combo29", 3, "thm:count29 (16.3 million triple boxes)",
        io -> Combo.run(cert("combo", "combo29_d8.npz"), 8, io; case=cert("combo", "cases", "case29_all.json"), mg2=5e-6, mg3=1e-6, mgm=1e-6, expect_bound=3.101462))
    add("profile_levels", 1, "profile-level bound max G = 3.2888648 at (22,0,3,0,0,0,0) and density < 0.631496 (anc/profile_levels.py; exact rationals on constants rounded up to 1e-13; Arb for the density)",
        io -> ProfileLevels.run(io; log26=joinpath(DATA, "logs", "density", "recheck_combo26_full.log")))
    tlog = tight === nothing ? joinpath(DATA, "logs", "density", "combo26_tight_check.log") : tight
    add("profile_levels_tight", 1, "prop:profile-levels, thm:density-cells: max G = 3.2826739 at (22,0,3,0,0,0,0) and density < 0.630996, with the pair and triple thresholds of the M = 26 run with pair margin 2e-7 (data/logs/density/combo26_tight_check.log; proved by that run and by the full-tier entry combo26_tight)",
        io -> ProfileLevels.run(io; tightlog=tlog))
    add("combo26_tight", 3, "the M = 26 typed certificate with pair margin 2e-7: the thresholds used by profile_levels_tight (7.5 million triple boxes)",
        io -> Combo.run(cert("combo", "combo26_level_d8.npz"), 8, io; case=cert("combo", "cases", "case26_all.json"), level="3.3352", mg2=2e-7, mg3=2e-6, mgm=2e-6, expect_bound=3.322615))
    if llm !== nothing
        add("explicit_eps0", 2, "thm:near-contact, thm:m24-robust (explicit eps_0)", io -> LLM24.run_eps0(llm, io))
        add("llm24_certificate", 2, "prop:verified, thm:twentyfour-points (steps 1, 2, 4, 6, 7 of the verification)",
            io -> LLM24.run_check(llm, io; expected_p2=cert("certificates", "llm24_p2.txt")))
    else
        add("explicit_eps0", 2, "thm:near-contact, thm:m24-robust (explicit eps_0)  -- needs --llm24 or --fetch-llm24", io -> nothing)
        add("llm24_certificate", 2, "prop:verified, thm:twentyfour-points  -- needs --llm24 or --fetch-llm24", io -> nothing)
    end
    E
end

function main(args)
    o = parse_args(args)
    llm = o["llm24"]
    if o["fetch"] !== nothing
        llm = LLM24.fetch(o["fetch"])
    end
    E = entries(llm, o["tight"])
    tier = TIERS[o["tier"]]
    sel = [e for e in E if e.tier <= tier && (o["only"] === nothing || occursin(o["only"], e.name))]
    if o["list"]
        for e in E; @printf("%-28s %-8s %s\n", e.name, ("quick", "standard", "full")[e.tier], e.supports); end
        return 0
    end
    mkpath(LOGS)
    rows = Tuple{String,String,Int,Int,Float64,String}[]
    for e in sel
        print(rpad(e.name, 28), " ... "); flush(stdout)
        t0 = time(); status = "PASS"; np = 0; nf = 0
        try
            R = open(joinpath(LOGS, e.name * ".log"), "w") do io
                e.f(io)
            end
            if R === nothing
                status = "SKIPPED"
            else
                np = count(c -> c.pass, R.checks); nf = count(c -> !c.pass, R.checks)
                status = nf == 0 ? "PASS" : "FAIL"
                if nf > 0
                    for c in R.checks; c.pass || println("\n    FAIL: ", c.name, "  got ", c.got, c.expected == "" ? "" : "; expected " * c.expected); end
                end
            end
        catch err
            status = "ERROR"
            open(joinpath(LOGS, e.name * ".log"), "a") do io; showerror(io, err, catch_backtrace()); end
        end
        dt = time() - t0
        @printf("%-7s %3d checks passed, %d failed  [%.0f s]\n", status, np, nf, dt)
        push!(rows, (e.name, status, np, nf, dt, e.supports))
    end
    open(joinpath(LOGS, o["only"] === nothing ? "SUMMARY.txt" : "SUMMARY_only.txt"), "w") do io
        for out in (stdout, io)
            println(out, "\n", rpad("check", 28), rpad("result", 9), rpad("passed", 8), rpad("failed", 8), rpad("seconds", 9), "supports")
            for (n, s, np, nf, dt, sup) in rows
                println(out, rpad(n, 28), rpad(s, 9), rpad(np, 8), rpad(nf, 8), rpad(round(Int, dt), 9), sup)
            end
            @printf(out, "total: %d run, %d passed, %d failed or in error, %d skipped; %.0f s\n", length(rows), count(r -> r[2] == "PASS", rows),
                    count(r -> r[2] in ("FAIL", "ERROR"), rows), count(r -> r[2] == "SKIPPED", rows), sum(r[5] for r in rows; init=0.0))
        end
    end
    any(r -> r[2] in ("FAIL", "ERROR"), rows) ? 1 : 0
end

exit(main(ARGS))
