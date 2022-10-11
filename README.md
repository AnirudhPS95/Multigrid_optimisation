# Multigrid_optimisation
2-D multigrid optimization subroutine in fortran

*******************************************************************
A simple poisson equation is solved using Gauss-Seidel scheme. 


The multigrid preconditioning is:

  1. One of the most powerful acceleration schemes
  2. It exploits the stability analysis of an iterative solver
  3. The longer wavelengths lead to slow convergence on fine grid
  4. Fine grid’s solution and coarse grid’s convergence was coupled
  5. This lead to development of multi-grid method
  
A V-goemtric multi grid was used to solve the following poission equation:

                  (𝜕^2 T)/(𝜕x^2 )+(𝜕^2 T)/(𝜕y^2 )=S

where, S (the source term) = −2(1−6𝑥^2 ) 𝑦^2 (1−𝑦^2 )+(1−6𝑦^2 ) 𝑥^2 (1−𝑥^2 )

The file with subroutine has two subroutines 
subroutine restriction(lvl,level,M,N,iter,dx,res,e,f) and
subroutine prolongation(lvl,level,M,N,iter,dx,e,f,er,T)






