#include <gss/formats/csv.hh>
#include <gss/homomorphism.hh>

#include <catch2/catch_test_macros.hpp>

#include <set>
#include <sstream>
#include <utility>

using namespace gss;

using std::chrono::operator""s;
using std::make_shared;
using std::make_unique;
using std::pair;
using std::set;
using std::string;
using std::stringstream;

TEST_CASE("homomorphism no edges")
{
    auto pattern = read_csv(stringstream{// clang-format off
R"(a,
b,
c,
)"}, "pattern"); // clang-format on

    auto target = read_csv(stringstream{// clang-format off
R"(1,
2,
3,
4,
)"}, "target"); // clang-format on

    HomomorphismParams params;
    params.timeout = make_shared<Timeout>(0s);
    params.restarts_schedule = make_unique<NoRestartsSchedule>();
    params.injectivity = Injectivity::NonInjective;

    SECTION("decide")
    {
        auto result = solve_homomorphism_problem(pattern, target, params);
        CHECK(result.mapping.size() == 3);
        CHECK(result.complete);
    }

    SECTION("count")
    {
        params.count_solutions = true;
        auto result = solve_homomorphism_problem(pattern, target, params);
        CHECK(result.solution_count == 4 * 4 * 4);
        CHECK(result.complete);
    }
}

TEST_CASE("homomorphism line")
{
    auto pattern = read_csv(stringstream{// clang-format off
R"(a,b
b,c
)"}, "pattern"); // clang-format on

    auto target = read_csv(stringstream{// clang-format off
R"(1,2
2,3
3,4
)"}, "target"); // clang-format on

    HomomorphismParams params;
    params.timeout = make_shared<Timeout>(0s);
    params.restarts_schedule = make_unique<NoRestartsSchedule>();
    params.injectivity = Injectivity::NonInjective;

    auto p = [&](const string & s) { return *pattern.vertex_from_name(s); };
    auto t = [&](const string & s) { return *target.vertex_from_name(s); };

    SECTION("count")
    {
        params.count_solutions = true;
        auto result = solve_homomorphism_problem(pattern, target, params);
        CHECK(result.solution_count == 10);
        CHECK(result.complete);
    }

    SECTION("enumerate")
    {
        set<VertexToVertexMapping> got;
        params.count_solutions = true;
        params.enumerate_callback = [&](const VertexToVertexMapping & m) -> bool {
            got.insert(m);
            return true;
        };
        auto result = solve_homomorphism_problem(pattern, target, params);
        CHECK(got == set<VertexToVertexMapping>{
                         {{p("a"), t("1")}, {p("b"), t("2")}, {p("c"), t("1")}}, //
                         {{p("a"), t("1")}, {p("b"), t("2")}, {p("c"), t("3")}}, //
                         {{p("a"), t("2")}, {p("b"), t("1")}, {p("c"), t("2")}}, //
                         {{p("a"), t("2")}, {p("b"), t("3")}, {p("c"), t("2")}}, //
                         {{p("a"), t("2")}, {p("b"), t("3")}, {p("c"), t("4")}}, //
                         {{p("a"), t("3")}, {p("b"), t("2")}, {p("c"), t("1")}}, //
                         {{p("a"), t("3")}, {p("b"), t("2")}, {p("c"), t("3")}}, //
                         {{p("a"), t("3")}, {p("b"), t("4")}, {p("c"), t("3")}}, //
                         {{p("a"), t("4")}, {p("b"), t("3")}, {p("c"), t("2")}}, //
                         {{p("a"), t("4")}, {p("b"), t("3")}, {p("c"), t("4")}} //
                     });
    }
}

TEST_CASE("homomorphism loop")
{
    auto pattern = read_csv(stringstream{// clang-format off
R"(a,a
a,b
)"}, "pattern"); // clang-format on

    auto target = read_csv(stringstream{// clang-format off
R"(1,2
2,2
2,3
3,4
)"}, "target"); // clang-format on

    HomomorphismParams params;
    params.timeout = make_shared<Timeout>(0s);
    params.restarts_schedule = make_unique<NoRestartsSchedule>();
    params.injectivity = Injectivity::NonInjective;

    auto p = [&](const string & s) { return *pattern.vertex_from_name(s); };
    auto t = [&](const string & s) { return *target.vertex_from_name(s); };

    SECTION("decide")
    {
        auto result = solve_homomorphism_problem(pattern, target, params);
        CHECK(result.solution_count == 1);
        CHECK(result.nodes == 0);
        CHECK(result.complete);
    }

    SECTION("count")
    {
        params.count_solutions = true;
        auto result = solve_homomorphism_problem(pattern, target, params);
        CHECK(result.solution_count == 3);
        CHECK(result.complete);
    }

    SECTION("enumerate")
    {
        set<VertexToVertexMapping> got;
        params.count_solutions = true;
        params.enumerate_callback = [&](const VertexToVertexMapping & m) -> bool {
            got.insert(m);
            return true;
        };
        auto result = solve_homomorphism_problem(pattern, target, params);
        CHECK(got == set<VertexToVertexMapping>{
                         {{p("a"), t("2")}, {p("b"), t("1")}}, //
                         {{p("a"), t("2")}, {p("b"), t("2")}}, //
                         {{p("a"), t("2")}, {p("b"), t("3")}} //
                     });
    }
}

TEST_CASE("induced mapping does not misuse the loop or clique shortcuts")
{
    // Regression: the non-injective target-loop shortcut and the clique-pattern shortcut both
    // ignored --induced on a loopy target, returning mappings that send a loopless pattern
    // vertex (or a pattern non-edge) onto a target self-loop, which an induced mapping forbids.

    SECTION("target-loop shortcut is not used for an induced mapping")
    {
        // Two isolated (loopless) pattern vertices; the only target vertex has a self-loop. The
        // non-injective loop shortcut would map both onto the loop, but no loopless pattern
        // vertex can map onto a looped target vertex under an induced mapping.
        auto pattern = read_csv(stringstream{"a,\nc,\n"}, "pattern");
        auto target = read_csv(stringstream{"1,1\n"}, "target");

        HomomorphismParams params;
        params.timeout = make_shared<Timeout>(0s);
        params.restarts_schedule = make_unique<NoRestartsSchedule>();
        params.injectivity = Injectivity::NonInjective;
        params.induced = true;

        auto decided = solve_homomorphism_problem(pattern, target, params);
        CHECK(decided.mapping.empty());
        CHECK(decided.complete);

        params.count_solutions = true;
        auto counted = solve_homomorphism_problem(pattern, target, params);
        CHECK(counted.solution_count == 0);
    }

    SECTION("clique shortcut is not used for an induced mapping into a loopy target")
    {
        // A K2 pattern is a simple (loopless) clique; the target's vertices all have self-loops.
        // The clique algorithm would find the K2, but its loopless endpoints cannot map onto
        // looped target vertices under an induced mapping.
        auto pattern = read_csv(stringstream{"a,b\n"}, "pattern");
        auto target = read_csv(stringstream{"1,2\n1,1\n2,2\n"}, "target");

        HomomorphismParams params;
        params.timeout = make_shared<Timeout>(0s);
        params.restarts_schedule = make_unique<NoRestartsSchedule>();
        params.clique_detection = true;
        params.induced = true;

        auto decided = solve_homomorphism_problem(pattern, target, params);
        CHECK(decided.mapping.empty());
        CHECK(decided.complete);

        params.count_solutions = true;
        auto counted = solve_homomorphism_problem(pattern, target, params);
        CHECK(counted.solution_count == 0);
    }

    SECTION("clique shortcut is still used for an induced mapping into a loopless target")
    {
        // The fix must not disable the optimisation in the common case: a K2 pattern into a
        // loopless K3 target has induced solutions and should still go through the clique
        // algorithm.
        auto pattern = read_csv(stringstream{"a,b\n"}, "pattern");
        auto target = read_csv(stringstream{"1,2\n1,3\n2,3\n"}, "target");

        HomomorphismParams params;
        params.timeout = make_shared<Timeout>(0s);
        params.restarts_schedule = make_unique<NoRestartsSchedule>();
        params.clique_detection = true;
        params.induced = true;

        auto decided = solve_homomorphism_problem(pattern, target, params);
        CHECK(! decided.mapping.empty());
        bool used_clique_solver = false;
        for (auto & stat : decided.extra_stats)
            if (stat == "used_clique_solver = true")
                used_clique_solver = true;
        CHECK(used_clique_solver);
    }
}

TEST_CASE("homomorphism loop shrinking")
{
    auto pattern = read_csv(stringstream{// clang-format off
R"(a,a
a,b
b,c
)"}, "pattern"); // clang-format on

    auto target = read_csv(stringstream{// clang-format off
R"(1,1
1,2
)"}, "target"); // clang-format on

    HomomorphismParams params;
    params.timeout = make_shared<Timeout>(0s);
    params.restarts_schedule = make_unique<NoRestartsSchedule>();
    params.injectivity = Injectivity::NonInjective;

    auto p = [&](const string & s) { return *pattern.vertex_from_name(s); };
    auto t = [&](const string & s) { return *target.vertex_from_name(s); };

    SECTION("decide")
    {
        auto result = solve_homomorphism_problem(pattern, target, params);
        CHECK(result.solution_count == 1);
        CHECK(result.nodes == 0);
        CHECK(result.complete);
    }

    SECTION("count")
    {
        params.count_solutions = true;
        auto result = solve_homomorphism_problem(pattern, target, params);
        CHECK(result.solution_count == 3);
        CHECK(result.complete);
    }

    SECTION("enumerate")
    {
        set<VertexToVertexMapping> got;
        params.count_solutions = true;
        params.enumerate_callback = [&](const VertexToVertexMapping & m) -> bool {
            got.insert(m);
            return true;
        };
        auto result = solve_homomorphism_problem(pattern, target, params);
        CHECK(got == set<VertexToVertexMapping>{
                         {{p("a"), t("1")}, {p("b"), t("1")}, {p("c"), t("1")}}, //
                         {{p("a"), t("1")}, {p("b"), t("1")}, {p("c"), t("2")}}, //
                         {{p("a"), t("1")}, {p("b"), t("2")}, {p("c"), t("1")}}, //
                     });
    }
}

TEST_CASE("minimum-cost homomorphism objective")
{
    SECTION("finds the lowest-cost structural mapping and bypasses the clique shortcut")
    {
        auto pattern = read_csv(stringstream{"a,b\n"}, "pattern");
        auto target = read_csv(stringstream{"1,2\n3,4\n"}, "target");

        auto t = [&](const string & s) { return unsigned(*target.vertex_from_name(s)); };

        HomomorphismParams params;
        params.timeout = make_shared<Timeout>(0s);
        params.restarts_schedule = make_unique<NoRestartsSchedule>();

        HomomorphismObjective objective;
        objective.pair_cost = [&](unsigned p, unsigned q, unsigned x, unsigned y) {
            if (p == q)
                return 0.0;
            const bool expensive_edge =
                (x == t("1") && y == t("2")) || (x == t("2") && y == t("1"));
            const bool cheap_edge =
                (x == t("3") && y == t("4")) || (x == t("4") && y == t("3"));
            if (expensive_edge)
                return 10.0;
            if (cheap_edge)
                return 1.0;
            return 1000.0;
        };
        params.objective = objective;

        auto result = solve_homomorphism_problem(pattern, target, params);
        REQUIRE(result.objective_value);
        CHECK(*result.objective_value == 1.0);
        CHECK(result.mapping.size() == 2);
        CHECK(result.complete);

        set<int> mapped_targets;
        for (auto & [_, target_vertex] : result.mapping)
            mapped_targets.insert(target_vertex);
        CHECK(mapped_targets == set<int>{int(t("3")), int(t("4"))});

        for (auto & stat : result.extra_stats)
            CHECK(stat != "used_clique_solver = true");
    }

    SECTION("uses the initial upper bound strictly and scores propagation-forced assignments")
    {
        auto pattern = read_csv(stringstream{"a,b\n"}, "pattern");
        auto target = read_csv(stringstream{"1,2\n3,4\n"}, "target");

        auto t = [&](const string & s) { return unsigned(*target.vertex_from_name(s)); };

        HomomorphismParams params;
        params.timeout = make_shared<Timeout>(0s);
        params.restarts_schedule = make_unique<NoRestartsSchedule>();

        HomomorphismObjective objective;
        objective.pair_cost = [&](unsigned p, unsigned q, unsigned x, unsigned y) {
            if (p == q)
                return 0.0;
            const bool expensive_edge =
                (x == t("1") && y == t("2")) || (x == t("2") && y == t("1"));
            const bool cheap_edge =
                (x == t("3") && y == t("4")) || (x == t("4") && y == t("3"));
            if (expensive_edge)
                return 10.0;
            if (cheap_edge)
                return 1.0;
            return 1000.0;
        };
        objective.initial_upper_bound = 1.0;
        params.objective = objective;

        auto result = solve_homomorphism_problem(pattern, target, params);
        CHECK(result.mapping.empty());
        CHECK_FALSE(result.objective_value);
        CHECK(result.complete);
    }

    SECTION("rejects semantically invalid assignments")
    {
        auto pattern = read_csv(stringstream{"a,b\n"}, "pattern");
        auto target = read_csv(stringstream{"1,2\n3,4\n"}, "target");

        auto t = [&](const string & s) { return unsigned(*target.vertex_from_name(s)); };

        HomomorphismParams params;
        params.timeout = make_shared<Timeout>(0s);
        params.restarts_schedule = make_unique<NoRestartsSchedule>();

        HomomorphismObjective objective;
        objective.vertex_allowed = [&](unsigned, unsigned tv) {
            return tv != t("1") && tv != t("2");
        };
        objective.pair_cost = [&](unsigned p, unsigned q, unsigned, unsigned) {
            return p == q ? 0.0 : 1.0;
        };
        params.objective = objective;

        auto result = solve_homomorphism_problem(pattern, target, params);
        REQUIRE(result.objective_value);
        set<int> mapped_targets;
        for (auto & [_, target_vertex] : result.mapping)
            mapped_targets.insert(target_vertex);
        CHECK(mapped_targets == set<int>{int(t("3")), int(t("4"))});
    }

    SECTION("respects objective call limits before the first solution")
    {
        auto pattern = read_csv(stringstream{"a,b\n"}, "pattern");
        auto target = read_csv(stringstream{"1,2\n3,4\n"}, "target");

        HomomorphismParams params;
        params.timeout = make_shared<Timeout>(0s);
        params.restarts_schedule = make_unique<NoRestartsSchedule>();

        HomomorphismObjective objective;
        objective.first_solution_call_limit = 1;
        objective.after_solution_call_limit = 1;
        objective.pair_cost = [](unsigned p, unsigned q, unsigned, unsigned) {
            return p == q ? 0.0 : 1.0;
        };
        params.objective = objective;

        auto result = solve_homomorphism_problem(pattern, target, params);
        CHECK(result.mapping.empty());
        CHECK_FALSE(result.complete);
    }

    SECTION("limits the number of accepted improving complete mappings")
    {
        auto pattern = read_csv(stringstream{"a,b\n"}, "pattern");
        auto target = read_csv(stringstream{"1,2\n3,4\n"}, "target");

        auto t = [&](const string & s) { return unsigned(*target.vertex_from_name(s)); };

        HomomorphismParams params;
        params.timeout = make_shared<Timeout>(0s);
        params.restarts_schedule = make_unique<NoRestartsSchedule>();

        params.value_ordering_heuristic = ValueOrdering::None;
        HomomorphismObjective objective;
        objective.max_solution_trials = 1;
        objective.pair_cost = [&](unsigned p, unsigned q, unsigned x, unsigned y) {
            if (p == q)
                return 0.0;
            const bool first_edge =
                (x == t("1") && y == t("2")) || (x == t("2") && y == t("1"));
            return first_edge ? 10.0 : 1.0;
        };
        params.objective = objective;

        auto result = solve_homomorphism_problem(pattern, target, params);
        REQUIRE(result.objective_value);
        CHECK(result.mapping.size() == 2);
        CHECK_FALSE(result.complete);
    }

    SECTION("includes vertex costs")
    {
        auto pattern = read_csv(stringstream{"a,\nb,\n"}, "pattern");
        auto target = read_csv(stringstream{"1,\n2,\n3,\n"}, "target");

        auto p = [&](const string & s) { return unsigned(*pattern.vertex_from_name(s)); };
        auto t = [&](const string & s) { return unsigned(*target.vertex_from_name(s)); };

        HomomorphismParams params;
        params.timeout = make_shared<Timeout>(0s);
        params.restarts_schedule = make_unique<NoRestartsSchedule>();

        HomomorphismObjective objective;
        objective.vertex_cost = [&](unsigned pv, unsigned tv) {
            if (pv == p("a"))
                return tv == t("1") ? 0.25 : 5.0;
            if (pv == p("b"))
                return tv == t("2") ? 0.5 : 5.0;
            return 0.0;
        };
        objective.order_values_by_cost = true;
        params.objective = objective;

        auto result = solve_homomorphism_problem(pattern, target, params);
        REQUIRE(result.objective_value);
        CHECK(*result.objective_value == 0.75);
        CHECK(result.mapping.at(p("a")) == t("1"));
        CHECK(result.mapping.at(p("b")) == t("2"));
        CHECK(result.complete);
    }
}
