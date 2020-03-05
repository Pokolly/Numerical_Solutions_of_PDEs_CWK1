%The basis of the FTCS method:
%Find The number of time steps, number of space steps, delta t and delta x
%Then use the boundary conditions to find U(0,0) & U(0,N)
%Then use the intial conditino to find U(1,j)
%Then repeat these steps for all the other time steps

boundarycondition_x_1 = -1; %set up our boundary conditions for x values
boundarycondition_x_2 = 1;
number_of_timesteps=50 ;%time steps 
N=20; %to account for matlab it's N-1 really
deltax = 1/20; %set up the spacings for delta x and delta t
deltat = 0.001;
r = deltat / ((deltax)^2); %pre-calculate r
M = ones(number_of_timesteps,N+1); %set up a 50x20 array of ones to store our solutions
for i=1:number_of_timesteps %places our 1st initial condition for x in the first column
    M(i,1) = boundarycondition_x_1 ; 
end
    
for i=1:number_of_timesteps %places our final boundary conditions for x in the final column
    M(i,N+1) = boundarycondition_x_2 ;
end

%calculate the first row of our matrix by using the formula given
for i = 1:N-1 
    M(1,i+1) = 2*deltax*i-1 ; 
end

%now for the piece de resistance, we then use our formula
%which is the PDE for iterating to calculate the rest of the terms in
%the matrix. 
for k=1:number_of_timesteps-1 %TIME STEP: we start at 1 and end with the number of timestamps -1
     for i=2:N  %SPACE STEP: noting that we start at 2 as we already have the first one complete
           M(k+1,i) = M(k,i)+ r*(M(k,i-1)-2*M(k,i)+M(k,i+1))+deltat*50*M(k,i)*(1-(M(k,i)^2));
     end
end
%Now we have our M matrix with all of our nice approximations for U in it.
%Please note that the columns are the x spacings and the rows are the time
%spacings.
surf(M) %creates a 3D surface of the Matrix M
%NOTE: my labelled version of the surface can be found uploaded to SL



