# Multigrid_optimisation
2-D multigrid optimization subroutine in fortran

*******************************************************************
A simple poisson equation is solved using Gauss-Seidel scheme. 


The multigrid preconditioning is:

  1. One of the most powerful acceleration schemes
  2. It exploits the stability analysis of an iterative solver
  3. The longer wavelengths lead to slow convergence on fine grid
  4. Fine gridβs solution and coarse gridβs convergence was coupled
  5. This lead to development of multi-grid method
  
A V-goemtric multi grid was used to solve the following poission equation:

                  (π^2 T)/(πx^2 )+(π^2 T)/(πy^2 )=S

where, S (the source term) = β2(1β6π₯^2 ) π¦^2 (1βπ¦^2 )+(1β6π¦^2 ) π₯^2 (1βπ₯^2 )

The file with subroutine has two subroutines 
subroutine restriction(lvl,level,M,N,iter,dx,res,e,f) and
subroutine prolongation(lvl,level,M,N,iter,dx,e,f,er,T)






