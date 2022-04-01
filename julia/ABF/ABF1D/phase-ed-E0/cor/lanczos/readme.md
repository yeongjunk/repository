2022/3/21
1D nu = 2 weak phase disordered ABF, scale free model.
When the system size is odd, it has strict E = 0 state

It finds the E = 0 states using sparse exact diagonalization and calculate IPR of that state.
Don't let E_c too small because there will be high proabbilty to fail if the L is large
It is inherently unstable because E_c may hit the actual eigenvalue.
This has to be fixed soon.
Also, nev should be larger than some value, which depends on E_c. Try finding it with experiment with small L and scale it with L


