use system
use ellipsoid

implicit none
integer j, i

print*, 'Program Crystal -- Hamaker'
print*, 'Calculates the Particle-Particle Interaction Energy'
print*, 'GIT Version: ', _VERSION
call readinput
call inittransf ! Create transformation matrixes
call calcHam ! calculates Hamaker interaction between particles
end


