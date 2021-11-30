import argparse
from project import simplex_step, simplex_init, simplex_method
from test_scripts import test1, test2, test3, test4, test5, test6, test7, test8, test9, test10, test11

parser = argparse.ArgumentParser()
parser.add_argument("--all", help="Run all test cases", action="store_true")
parser.add_argument("--test", metavar='', help="Run one specified test case", type=int)
parser.add_argument("--function", metavar='', help="Run the appropriate tests for a specified function", type=str)

# Organize test cases
step_cases = [test1, test2, test3, test4, test5]
init_cases = [test6, test7, test8]
method_cases = [test9, test10, test11]


def validateCase(_exec):
    if _exec:
        print("Test successful!")


if __name__ == "__main__":
    args = parser.parse_args()
    if args.all:
        # Run all cases >>
        for case in step_cases:
            validateCase(case.run(simplex_step))
        for case in init_cases:
            validateCase(case.run(simplex_init))
        for case in method_cases:
            validateCase(case.run(simplex_method))
    elif args.test:
        assert (args.test >= 1 and args.test <= 9), "Invalid test number!"
        if args.test in range(1, 6):
            validateCase(step_cases[args.test - 1].run(simplex_step))
        elif args.test in range(6, 9):
            validateCase(init_cases[args.test - 6].run(simplex_init))
        else:
            validateCase(method_cases[args.test - 9].run(simplex_method))
    elif args.function:
        assert (args.function in ["simplex_step", "simplex_init", "simplex_method"]), "Invalid function!"
        if args.function == "simplex_step":
            for case in step_cases:
                validateCase(case.run(simplex_step))
        elif args.function == "simplex_init":
            for case in init_cases:
                validateCase(case.run(simplex_init))
        else:
            for case in method_cases:
                validateCase(case.run(simplex_method))
    else:
        assert (False), "GROOT doesn't have any instruction!"
