![image](https://cdn.dribbble.com/users/4353062/screenshots/8179725/media/b823df3cf607d39fe77d44b0ecf56105.jpg?compress=1&resize=400x300)

![image](https://img.shields.io/badge/Version-1.0-green)

*Making Simplex Project Testing - Clean and Simple*

- **What does this repo do?** 
  - It organizes the python stack for the coding project
  

- **What do I need to do in the project?**
  - The project aim is to complete 3 functions using the revised simplex method
  - These are `simplex_step`, `simplex_init`, and `simplex_method`
  - The details of what these functions do and what are the inputs and outputs are available in project.pdf file provided.
  - The instructors will run the functions for different test cases (1-11)
  

- **Okay, then how is this repository related to project?**
  - Find `project.py` file. This is the only file where you make changes.
  - You'll find those 3 functions in this file, just add the appropriate code in there
  - Once done, you can test your implementation using `run_test.py` and verify if your implementation is correct.
  - Every test file provided by the instructor tests different aspects of the simplex implementation, so to make things simple, the `run_test.py` by itself will find the appropriate test case - run your function and then let you know the result (success/fail)
  - There's more, groot also checks if your function return types are valid, so that you know that you have the right construction
  - Beyond that, groot is very flexible in terms of testing, you can test by test_number, function or all cases in one go. 
  

- **"Uhhh.. that seems to complicated, tell me bottom line?"**
  - Just go to project.py file and write your code.
  - Say, you've finished `simplex_step` and want to see if this is correct.
    - On your console run - `python3 run_test.py --function=simplex_step`
  - OR say you just want to run test3 for debugging sake
    - On your console run - `python3 run_test.py --test=3`
    - OR may be you're done with all 3 functions and run all the test cases
      - On your console run - `python3 run_test.py --all`
  

  - **Nice! What do I get after running this?**
    - Groot tells you if your output were correct/incorrect and if latter then what part of the output is incorrect. Like this:
    ```
    $python run_test.py --function="simplex_step"
    
    ************************************
    TEST-1
    ************************************
    INCORRECT STEP!
    iN incorrect!
    iB incorrect!
    
    >> TEST - 2 was successful! 
    >> TEST - 3 was successful!
    >> TEST - 4 was successful!
    ************************************
    TEST-5
    ************************************
    INCORRECT istatus!
    ```

