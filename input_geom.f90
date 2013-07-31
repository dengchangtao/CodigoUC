!Rutina que lee la batimetria y la guarda en x_global, y_global, z_global

SUBROUTINE input_geom

USE global_variables
USE geometries

implicit none

integer::i,j

allocate (x_global(Nbx,Nby),y_global(Nbx,Nby),z_global(Nbx,Nby))

print*, Nbx, Nby
open	(unit=2, file ='gridX.dat')!, form='unformatted')
read	(2,*) ((x_global(i,j),i=1,Nbx),j=1,Nby)
close(unit=2)

open	(unit=3, file ='gridY.dat')!, form='unformatted')
read	(3,*) ((y_global(i,j),i=1,Nbx),j=1,Nby)
close(unit=3)

open	(unit=4,file='gridZ.dat')!, form='unformatted')
read	(4,*) ((z_global(i,j),i=1,Nbx),j=1,Nby)
close(unit=4)



END SUBROUTINE input_geom
