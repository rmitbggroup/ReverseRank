Install

The Installation of GLPK is required. 
Steps:
1. Select the packages regarding your operating system. (win/linux)
2. Add source code of glpk in to cmake including directory.
    include_directories($address_of_source$)
3. Target the link library
    TARGET_LINK_LIBRARIES($project$ $address_of_lib$)
Please take the CMakeLists.txt as an example.

Running

1. Setting of the parameters
    a. N, number of dimensions
    b. M, number of the items
    c. path, the root directory for files to load in
    d. tp, the location of transition graph.
    e. awp, the location of random initialized vector
    f. itemp, the location of all the items
    g. wgp, the location of files with top-100 users
2. Run the tests, start from VaryTest(). 
3. Test different methods from the input parameters of VaryTest().
    a. set the method. 1 for lp, 2 for itrssp, 3 for hybrid, 4 for greedy, 5 for itrlp.
    b. set the rounds, which lead to how many times this method going to operate.
    c. set the ue = 1, if you want to record the results.
4. Note that, the results are located in the same address with path.


    