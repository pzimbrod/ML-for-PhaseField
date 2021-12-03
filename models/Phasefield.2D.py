from math import pi,cos,tanh,sqrt,pi

def set_left_boundary(field,mgrid,ngrid):
        for k in range(ngrid):
                field[0][k] = field[1][k]
#******************************************************************
def set_right_boundary(field,mgrid,ngrid):
        for k in range(ngrid):
                field[mgrid-1][k] = field[mgrid-2][k]
#******************************************************************
def set_top_boundary(field,mgrid,ngrid):
        for i in range(mgrid):
                field[i][ngrid-1] = field[i][ngrid-2]
#******************************************************************
def set_bottom_boundary(field,mgrid,ngrid):
        for i in range(mgrid):
                field[i][0] = field[i][1]
#******************************************************************
def gprim(phi):
        return 2*phi*(1.0-phi)*(1.0-2.0*phi)
#******************************************************************
def hprim(phi):
        return 6.0*phi*(1.0-phi)
#******************************************************************
def initialize_system(phi,delta_x,radius,mgrid,ngrid):
        for i in range(mgrid):
                for k in range(ngrid):
                        if (sqrt((i-0.0*ngrid/2.)**2 + (k-0.0*mgrid/2.)**2)*delta_x <= radius):
                                phi[i][k] = 1.0
                        else:
                                phi[i][k] = 0.0
#******************************************************************
def solve_phasefield(phi,dphi,delta_x,delta_t,xi,M,gamma_s,mu0,mgrid,ngrid):
        for i in range(1,mgrid-1):
                for k in range(1,ngrid-1):# for every inner point
                        laplace=(phi[i+1][k] + phi[i-1][k] + phi[i][k+1] + phi[i][k-1] - 4.0*phi[i][k]) / (delta_x*delta_x)
                        dphi[i][k] =  laplace - 2.0*gprim(phi[i][k])/(xi*xi) + mu0/(3*gamma_s*xi)*hprim(phi[i][k])
                        phi[i][k] += M * delta_t * dphi[i][k]
                                
        set_left_boundary(phi,mgrid,ngrid)
        set_right_boundary(phi,mgrid,ngrid)
        set_top_boundary(phi,mgrid,ngrid)
        set_bottom_boundary(phi,mgrid,ngrid)

#******************************************************************
        
def save_field(field,filename):
        outputfile = open(filename,'w')
        for i in range(len(field)):
                for k in range(len(field[i])):
                        outputfile.write(str(i) + " " + str(k)  + " " + str(field[i][k]))
                outputfile.write(' ')
        outputfile.close()

#******************************************************************

def integral(field, delta_x):
        integral = 0
        for i in range(1,len(field)):
                for k in range(1,len(field[i])):
                        integral += field[i][k] * delta_x**2
        return integral

#******************************************************************


# Initialize parameters (time step, grid spacing, etc.).

mgrid = 55 # number of grid points in x direction, including boundary points
ngrid = 50 # number of grid points in y direction, including boundary points

delta_t = 1.e-1                 # [t] timestep 
nloop = 1001             # number of iterations 
  
delta_x = 1.0                # [l] physical lattice parameter 
  
gams = 50.0                 # [m/t^2] interface energy 
xi   = 2.5                  # [l] phasefield interface width 
M    = 1.0                  # [l^2/t] kinetic coefficient 
mu0  = 0.0			        # mu0 = m * (T-TM) 
initial_radius = 20

outstep_profile=100 # fields written out every outstep steps
outstep_radius=10
time = 0.0
output_radius = open("radius.dat","w")

#Initialize fields
phi = [[0.0 for k in range(ngrid)] for i in range(mgrid)]
dphi = [[0.0 for k in range(ngrid)] for i in range(mgrid)]
  
initialize_system(phi,delta_x,initial_radius,mgrid,ngrid)

print('\n *** Start the time loop *** \n')
	
for Step in range(nloop):
        
        solve_phasefield(phi,dphi,delta_x,delta_t,xi,M,gams,mu0,mgrid,ngrid)
       
        time += delta_t
        
# Output
        if Step%outstep_radius == 0: # Every outstep_radius steps
        
                #write out approximate radius of nucleus
                area = integral(phi,delta_x)
                radius = sqrt(4.0*area/pi)
                output_radius.write(str(time) + '   ' + str(radius) + ' \n')
                print('time = ', time, 'radius = ', radius)
		
        if Step%outstep_profile == 0: # Every outstep_profile steps
        
                filenamePhi = str(Step/outstep_profile + 1)+".phi.dat"
                
                save_field(phi,filenamePhi)
		
print('\n *** Program finished *** \n')
