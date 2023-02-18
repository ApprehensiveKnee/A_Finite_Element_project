%% Define the x, y, z coordinates of the points
clc;
clear all;
close all;

x = [0 0.1 0 0.1 0 0.1 0 0.1 0 0.1 0 0.1 0 0.1 0 0.1 0 0.1 0 0.1 0 0.1 0 0.1 0 0.1 0 0.1 0 0.1 0 0.1 0 0.1 0 0.1 0 0.1 0 0.1 0.1 0.2 0.1 0.2 0.1 0.2 0.1 0.2 0.1 0.2 0.1 0.2 0.1 0.2 0.1 0.2 0.1 0.2 0.1 0.2 0.1 0.2 0.1 0.2 0.1 0.2 0.1 0.2 0.1 0.2 0.1 0.2 0.1 0.2 0.1 0.2 0.1 0.2 0.1 0.2 0.2 0.3 0.2 0.3 0.2 0.3 0.2 0.3 0.2 0.3 0.2 0.3 0.2 0.3 0.2 0.3 0.2 0.3 0.2 0.3 0.2 0.3 0.2 0.3 0.2 0.3 0.2 0.3 0.2 0.3 0.2 0.3 0.2 0.3 0.2 0.3 0.2 0.3 0.2 0.3 0.3 0.4 0.3 0.4 0.3 0.4 0.3 0.4 0.3 0.4 0.3 0.4 0.3 0.4 0.3 0.4 0.3 0.4 0.3 0.4 0.3 0.4 0.3 0.4 0.3 0.4 0.3 0.4 0.3 0.4 0.3 0.4 0.3 0.4 0.3 0.4 0.3 0.4 0.3 0.4 0.4 0.5 0.4 0.5 0.4 0.5 0.4 0.5 0.4 0.5 0.4 0.5 0.4 0.5 0.4 0.5 0.4 0.5 0.4 0.5 0.4 0.5 0.4 0.5 0.4 0.5 0.4 0.5 0.4 0.5 0.4 0.5 0.4 0.5 0.4 0.5 0.4 0.5 0.4 0.5 0.5 0.6 0.5 0.6 0.5 0.6 0.5 0.6 0.5 0.6 0.5 0.6 0.5 0.6 0.5 0.6 0.5 0.6 0.5 0.6 0.5 0.6 0.5 0.6 0.5 0.6 0.5 0.6 0.5 0.6 0.5 0.6 0.5 0.6 0.5 0.6 0.5 0.6 0.5 0.6 0.6 0.7 0.6 0.7 0.6 0.7 0.6 0.7 0.6 0.7 0.6 0.7 0.6 0.7 0.6 0.7 0.6 0.7 0.6 0.7 0.6 0.7 0.6 0.7 0.6 0.7 0.6 0.7 0.6 0.7 0.6 0.7 0.6 0.7 0.6 0.7 0.6 0.7 0.6 0.7 0.7 0.8 0.7 0.8 0.7 0.8 0.7 0.8 0.7 0.8 0.7 0.8 0.7 0.8 0.7 0.8 0.7 0.8 0.7 0.8 0.7 0.8 0.7 0.8 0.7 0.8 0.7 0.8 0.7 0.8 0.7 0.8 0.7 0.8 0.7 0.8 0.7 0.8 0.7 0.8 0.8 0.9 0.8 0.9 0.8 0.9 0.8 0.9 0.8 0.9 0.8 0.9 0.8 0.9 0.8 0.9 0.8 0.9 0.8 0.9 0.8 0.9 0.8 0.9 0.8 0.9 0.8 0.9 0.8 0.9 0.8 0.9 0.8 0.9 0.8 0.9 0.8 0.9 0.8 0.9 0.9 1 0.9 1 0.9 1 0.9 1 0.9 1 0.9 1 0.9 1 0.9 1 0.9 1 0.9 1 0.9 1 0.9 1 0.9 1 0.9 1 0.9 1 0.9 1 0.9 1 0.9 1 0.9 1 0.9 1 ]

y = [0 0 0.1 0.1 0.1 0.1 0.2 0.2 0.2 0.2 0.3 0.3 0.3 0.3 0.4 0.4 0.4 0.4 0.5 0.5 0.5 0.5 0.6 0.6 0.6 0.6 0.7 0.7 0.7 0.7 0.8 0.8 0.8 0.8 0.9 0.9 0.9 0.9 1 1 0 0 0.1 0.1 0.1 0.1 0.2 0.2 0.2 0.2 0.3 0.3 0.3 0.3 0.4 0.4 0.4 0.4 0.5 0.5 0.5 0.5 0.6 0.6 0.6 0.6 0.7 0.7 0.7 0.7 0.8 0.8 0.8 0.8 0.9 0.9 0.9 0.9 1 1 0 0 0.1 0.1 0.1 0.1 0.2 0.2 0.2 0.2 0.3 0.3 0.3 0.3 0.4 0.4 0.4 0.4 0.5 0.5 0.5 0.5 0.6 0.6 0.6 0.6 0.7 0.7 0.7 0.7 0.8 0.8 0.8 0.8 0.9 0.9 0.9 0.9 1 1 0 0 0.1 0.1 0.1 0.1 0.2 0.2 0.2 0.2 0.3 0.3 0.3 0.3 0.4 0.4 0.4 0.4 0.5 0.5 0.5 0.5 0.6 0.6 0.6 0.6 0.7 0.7 0.7 0.7 0.8 0.8 0.8 0.8 0.9 0.9 0.9 0.9 1 1 0 0 0.1 0.1 0.1 0.1 0.2 0.2 0.2 0.2 0.3 0.3 0.3 0.3 0.4 0.4 0.4 0.4 0.5 0.5 0.5 0.5 0.6 0.6 0.6 0.6 0.7 0.7 0.7 0.7 0.8 0.8 0.8 0.8 0.9 0.9 0.9 0.9 1 1 0 0 0.1 0.1 0.1 0.1 0.2 0.2 0.2 0.2 0.3 0.3 0.3 0.3 0.4 0.4 0.4 0.4 0.5 0.5 0.5 0.5 0.6 0.6 0.6 0.6 0.7 0.7 0.7 0.7 0.8 0.8 0.8 0.8 0.9 0.9 0.9 0.9 1 1 0 0 0.1 0.1 0.1 0.1 0.2 0.2 0.2 0.2 0.3 0.3 0.3 0.3 0.4 0.4 0.4 0.4 0.5 0.5 0.5 0.5 0.6 0.6 0.6 0.6 0.7 0.7 0.7 0.7 0.8 0.8 0.8 0.8 0.9 0.9 0.9 0.9 1 1 0 0 0.1 0.1 0.1 0.1 0.2 0.2 0.2 0.2 0.3 0.3 0.3 0.3 0.4 0.4 0.4 0.4 0.5 0.5 0.5 0.5 0.6 0.6 0.6 0.6 0.7 0.7 0.7 0.7 0.8 0.8 0.8 0.8 0.9 0.9 0.9 0.9 1 1 0 0 0.1 0.1 0.1 0.1 0.2 0.2 0.2 0.2 0.3 0.3 0.3 0.3 0.4 0.4 0.4 0.4 0.5 0.5 0.5 0.5 0.6 0.6 0.6 0.6 0.7 0.7 0.7 0.7 0.8 0.8 0.8 0.8 0.9 0.9 0.9 0.9 1 1 0 0 0.1 0.1 0.1 0.1 0.2 0.2 0.2 0.2 0.3 0.3 0.3 0.3 0.4 0.4 0.4 0.4 0.5 0.5 0.5 0.5 0.6 0.6 0.6 0.6 0.7 0.7 0.7 0.7 0.8 0.8 0.8 0.8 0.9 0.9 0.9 0.9 1 1 ]
z = [0 0 0 0.559017 0 0.559017 0 0.345492 0 0.345492 0 -0.345492 0 -0.345492 0 -0.559017 0 -0.559017 0 0 0 0 0 0.559017 0 0.559017 0 0.345492 0 0.345492 0 -0.345492 0 -0.345492 0 -0.559017 0 -0.559017 0 0 0 0 0.559017 0.904508 0.559017 0.904508 0.345492 0.559017 0.345492 0.559017 -0.345492 -0.559017 -0.345492 -0.559017 -0.559017 -0.904508 -0.559017 -0.904508 0 0 0 0 0.559017 0.904508 0.559017 0.904508 0.345492 0.559017 0.345492 0.559017 -0.345492 -0.559017 -0.345492 -0.559017 -0.559017 -0.904508 -0.559017 -0.904508 0 0 0 0 0.904508 0.904508 0.904508 0.904508 0.559017 0.559017 0.559017 0.559017 -0.559017 -0.559017 -0.559017 -0.559017 -0.904508 -0.904508 -0.904508 -0.904508 0 0 0 0 0.904508 0.904508 0.904508 0.904508 0.559017 0.559017 0.559017 0.559017 -0.559017 -0.559017 -0.559017 -0.559017 -0.904508 -0.904508 -0.904508 -0.904508 0 0 0 0 0.904508 0.559017 0.904508 0.559017 0.559017 0.345492 0.559017 0.345492 -0.559017 -0.345492 -0.559017 -0.345492 -0.904508 -0.559017 -0.904508 -0.559017 0 0 0 0 0.904508 0.559017 0.904508 0.559017 0.559017 0.345492 0.559017 0.345492 -0.559017 -0.345492 -0.559017 -0.345492 -0.904508 -0.559017 -0.904508 -0.559017 0 0 0 0 0.559017 0 0.559017 0 0.345492 0 0.345492 0 -0.345492 0 -0.345492 0 -0.559017 0 -0.559017 0 0 0 0 0 0.559017 0 0.559017 0 0.345492 0 0.345492 0 -0.345492 0 -0.345492 0 -0.559017 0 -0.559017 0 0 0 0 0 0 -0.559017 0 -0.559017 0 -0.345492 0 -0.345492 0 0.345492 0 0.345492 0 0.559017 0 0.559017 0 0 0 0 0 -0.559017 0 -0.559017 0 -0.345492 0 -0.345492 0 0.345492 0 0.345492 0 0.559017 0 0.559017 0 0 0 0 -0.559017 -0.904508 -0.559017 -0.904508 -0.345492 -0.559017 -0.345492 -0.559017 0.345492 0.559017 0.345492 0.559017 0.559017 0.904508 0.559017 0.904508 0 0 0 0 -0.559017 -0.904508 -0.559017 -0.904508 -0.345492 -0.559017 -0.345492 -0.559017 0.345492 0.559017 0.345492 0.559017 0.559017 0.904508 0.559017 0.904508 0 0 0 0 -0.904508 -0.904508 -0.904508 -0.904508 -0.559017 -0.559017 -0.559017 -0.559017 0.559017 0.559017 0.559017 0.559017 0.904508 0.904508 0.904508 0.904508 0 0 0 0 -0.904508 -0.904508 -0.904508 -0.904508 -0.559017 -0.559017 -0.559017 -0.559017 0.559017 0.559017 0.559017 0.559017 0.904508 0.904508 0.904508 0.904508 0 0 0 0 -0.904508 -0.559017 -0.904508 -0.559017 -0.559017 -0.345492 -0.559017 -0.345492 0.559017 0.345492 0.559017 0.345492 0.904508 0.559017 0.904508 0.559017 0 0 0 0 -0.904508 -0.559017 -0.904508 -0.559017 -0.559017 -0.345492 -0.559017 -0.345492 0.559017 0.345492 0.559017 0.345492 0.904508 0.559017 0.904508 0.559017 0 0 0 0 -0.559017 0 -0.559017 0 -0.345492 0 -0.345492 0 0.345492 0 0.345492 0 0.559017 0 0.559017 0 0 0 0 0 -0.559017 0 -0.559017 0 -0.345492 0 -0.345492 0 0.345492 0 0.345492 0 0.559017 0 0.559017 0 0 0 ]
% Plot the values as a 3D scatter plot
scatter3(x, y, z)
xlabel('x')
ylabel('y')
zlabel('z')

%% E
clc;
clear all;
close all;
x = linspace(-1,1,2000);
z = [-0 -2.25067e+12 -4.4991e+12 -6.74527e+12 -8.98918e+12 -1.12309e+13 -1.34703e+13 -1.57074e+13 -1.79423e+13 -2.0175e+13 -2.24054e+13 -2.46336e+13 -2.68595e+13 -2.90831e+13 -3.13045e+13 -3.35237e+13 -3.57406e+13 -3.79552e+13 -4.01676e+13 -4.23777e+13 -4.45856e+13 -4.67913e+13 -4.89947e+13 -5.11958e+13 -5.33947e+13 -5.55913e+13 -5.77857e+13 -5.99778e+13 -6.21677e+13 -6.43553e+13 -6.65407e+13 -6.87238e+13 -7.09047e+13 -7.30833e+13 -7.52597e+13 -7.74338e+13 -7.96056e+13 -8.17752e+13 -8.39426e+13 -8.61077e+13 -8.82706e+13 -9.04312e+13 -9.25895e+13 -9.47456e+13 -9.68994e+13 -9.9051e+13 -1.012e+14 -1.03347e+14 -1.05492e+14 -1.07635e+14 -1.09775e+14 -1.11913e+14 -1.14049e+14 -1.16183e+14 -1.18314e+14 -1.20443e+14 -1.2257e+14 -1.24695e+14 -1.26817e+14 -1.28937e+14 -1.31055e+14 -1.3317e+14 -1.35284e+14 -1.37395e+14 -1.39504e+14 -1.4161e+14 -1.43714e+14 -1.45816e+14 -1.47916e+14 -1.50014e+14 -1.52109e+14 -1.54202e+14 -1.56293e+14 -1.58381e+14 -1.60468e+14 -1.62552e+14 -1.64634e+14 -1.66713e+14 -1.6879e+14 -1.70865e+14 -1.72938e+14 -1.75009e+14 -1.77077e+14 -1.79143e+14 -1.81207e+14 -1.83268e+14 -1.85328e+14 -1.87385e+14 -1.89439e+14 -1.91492e+14 -1.93542e+14 -1.9559e+14 -1.97636e+14 -1.99679e+14 -2.01721e+14 -2.0376e+14 -2.05796e+14 -2.07831e+14 -2.09863e+14 -2.11893e+14 -2.13921e+14 -2.15946e+14 -2.1797e+14 -2.19991e+14 -2.22009e+14 -2.24026e+14 -2.2604e+14 -2.28052e+14 -2.30062e+14 -2.32069e+14 -2.34075e+14 -2.36078e+14 -2.38078e+14 -2.40077e+14 -2.42073e+14 -2.44067e+14 -2.46059e+14 -2.48048e+14 -2.50035e+14 -2.5202e+14 -2.54003e+14 -2.55983e+14 -2.57962e+14 -2.59938e+14 -2.61911e+14 -2.63883e+14 -2.65852e+14 -2.67819e+14 -2.69784e+14 -2.71746e+14 -2.73706e+14 -2.75664e+14 -2.7762e+14 -2.79573e+14 -2.81525e+14 -2.83473e+14 -2.8542e+14 -2.87365e+14 -2.89307e+14 -2.91247e+14 -2.93184e+14 -2.9512e+14 -2.97053e+14 -2.98984e+14 -3.00913e+14 -3.02839e+14 -3.04763e+14 -3.06685e+14 -3.08605e+14 -3.10522e+14 -3.12437e+14 -3.1435e+14 -3.16261e+14 -3.18169e+14 -3.20075e+14 -3.21979e+14 -3.23881e+14 -3.2578e+14 -3.27677e+14 -3.29572e+14 -3.31465e+14 -3.33355e+14 -3.35243e+14 -3.37129e+14 -3.39013e+14 -3.40894e+14 -3.42773e+14 -3.4465e+14 -3.46525e+14 -3.48397e+14 -3.50267e+14 -3.52135e+14 -3.54001e+14 -3.55864e+14 -3.57725e+14 -3.59584e+14 -3.61441e+14 -3.63295e+14 -3.65147e+14 -3.66997e+14 -3.68845e+14 -3.7069e+14 -3.72533e+14 -3.74374e+14 -3.76213e+14 -3.78049e+14 -3.79883e+14 -3.81715e+14 -3.83545e+14 -3.85372e+14 -3.87197e+14 -3.8902e+14 -3.9084e+14 -3.92659e+14 -3.94475e+14 -3.96289e+14 -3.981e+14 -3.9991e+14 -4.01717e+14 -4.03521e+14 -4.05324e+14 -4.07124e+14 -4.08922e+14 -4.10718e+14 -4.12512e+14 -4.14303e+14 -4.16092e+14 -4.17879e+14 -4.19663e+14 -4.21446e+14 -4.23226e+14 -4.25004e+14 -4.26779e+14 -4.28552e+14 -4.30323e+14 -4.32092e+14 -4.33859e+14 -4.35623e+14 -4.37385e+14 -4.39145e+14 -4.40902e+14 -4.42658e+14 -4.44411e+14 -4.46161e+14 -4.4791e+14 -4.49656e+14 -4.514e+14 -4.53142e+14 -4.54882e+14 -4.56619e+14 -4.58354e+14 -4.60087e+14 -4.61817e+14 -4.63545e+14 -4.65271e+14 -4.66995e+14 -4.68717e+14 -4.70436e+14 -4.72153e+14 -4.73868e+14 -4.7558e+14 -4.7729e+14 -4.78998e+14 -4.80704e+14 -4.82408e+14 -4.84109e+14 -4.85808e+14 -4.87505e+14 -4.89199e+14 -4.90891e+14 -4.92581e+14 -4.94269e+14 -4.95954e+14 -4.97638e+14 -4.99319e+14 -5.00997e+14 -5.02674e+14 -5.04348e+14 -5.0602e+14 -5.0769e+14 -5.09357e+14 -5.11022e+14 -5.12685e+14 -5.14346e+14 -5.16004e+14 -5.17661e+14 -5.19315e+14 -5.20966e+14 -5.22616e+14 -5.24263e+14 -5.25908e+14 -5.27551e+14 -5.29191e+14 -5.30829e+14 -5.32465e+14 -5.34099e+14 -5.3573e+14 -5.37359e+14 -5.38986e+14 -5.40611e+14 -5.42233e+14 -5.43854e+14 -5.45471e+14 -5.47087e+14 -5.48701e+14 -5.50312e+14 -5.51921e+14 -5.53527e+14 -5.55132e+14 -5.56734e+14 -5.58334e+14 -5.59931e+14 -5.61527e+14 -5.6312e+14 -5.64711e+14 -5.663e+14 -5.67886e+14 -5.6947e+14 -5.71052e+14 -5.72632e+14 -5.74209e+14 -5.75784e+14 -5.77357e+14 -5.78928e+14 -5.80496e+14 -5.82062e+14 -5.83626e+14 -5.85188e+14 -5.86747e+14 -5.88304e+14 -5.89859e+14 -5.91412e+14 -5.92962e+14 -5.9451e+14 -5.96056e+14 -5.976e+14 -5.99141e+14 -6.0068e+14 -6.02217e+14 -6.03751e+14 -6.05284e+14 -6.06814e+14 -6.08342e+14 -6.09867e+14 -6.11391e+14 -6.12912e+14 -6.14431e+14 -6.15947e+14 -6.17462e+14 -6.18974e+14 -6.20483e+14 -6.21991e+14 -6.23496e+14 -6.24999e+14 -6.265e+14 -6.27999e+14 -6.29495e+14 -6.30989e+14 -6.32481e+14 -6.33971e+14 -6.35458e+14 -6.36943e+14 -6.38426e+14 -6.39906e+14 -6.41385e+14 -6.42861e+14 -6.44335e+14 -6.45806e+14 -6.47275e+14 -6.48742e+14 -6.50207e+14 -6.5167e+14 -6.5313e+14 -6.54588e+14 -6.56044e+14 -6.57497e+14 -6.58949e+14 -6.60398e+14 -6.61844e+14 -6.63289e+14 -6.64731e+14 -6.66171e+14 -6.67609e+14 -6.69045e+14 -6.70478e+14 -6.71909e+14 -6.73338e+14 -6.74764e+14 -6.76188e+14 -6.7761e+14 -6.7903e+14 -6.80448e+14 -6.81863e+14 -6.83276e+14 -6.84687e+14 -6.86095e+14 -6.87502e+14 -6.88906e+14 -6.90307e+14 -6.91707e+14 -6.93104e+14 -6.94499e+14 -6.95892e+14 -6.97282e+14 -6.9867e+14 -7.00056e+14 -7.0144e+14 -7.02822e+14 -7.04201e+14 -7.05578e+14 -7.06953e+14 -7.08325e+14 -7.09695e+14 -7.11063e+14 -7.12429e+14 -7.13792e+14 -7.15154e+14 -7.16513e+14 -7.17869e+14 -7.19224e+14 -7.20576e+14 -7.21926e+14 -7.23274e+14 -7.24619e+14 -7.25962e+14 -7.27303e+14 -7.28642e+14 -7.29978e+14 -7.31313e+14 -7.32644e+14 -7.33974e+14 -7.35302e+14 -7.36627e+14 -7.3795e+14 -7.3927e+14 -7.40589e+14 -7.41905e+14 -7.43219e+14 -7.44531e+14 -7.4584e+14 -7.47147e+14 -7.48452e+14 -7.49755e+14 -7.51055e+14 -7.52353e+14 -7.53649e+14 -7.54943e+14 -7.56234e+14 -7.57523e+14 -7.5881e+14 -7.60095e+14 -7.61377e+14 -7.62658e+14 -7.63935e+14 -7.65211e+14 -7.66485e+14 -7.67756e+14 -7.69025e+14 -7.70291e+14 -7.71556e+14 -7.72818e+14 -7.74078e+14 -7.75335e+14 -7.76591e+14 -7.77844e+14 -7.79095e+14 -7.80343e+14 -7.8159e+14 -7.82834e+14 -7.84076e+14 -7.85315e+14 -7.86553e+14 -7.87788e+14 -7.89021e+14 -7.90251e+14 -7.91479e+14 -7.92706e+14 -7.93929e+14 -7.95151e+14 -7.9637e+14 -7.97587e+14 -7.98802e+14 -8.00015e+14 -8.01225e+14 -8.02433e+14 -8.03639e+14 -8.04843e+14 -8.06044e+14 -8.07243e+14 -8.0844e+14 -8.09635e+14 -8.10827e+14 -8.12017e+14 -8.13205e+14 -8.1439e+14 -8.15574e+14 -8.16755e+14 -8.17934e+14 -8.1911e+14 -8.20285e+14 -8.21457e+14 -8.22626e+14 -8.23794e+14 -8.24959e+14 -8.26122e+14 -8.27283e+14 -8.28442e+14 -8.29598e+14 -8.30752e+14 -8.31904e+14 -8.33053e+14 -8.34201e+14 -8.35346e+14 -8.36488e+14 -8.37629e+14 -8.38767e+14 -8.39903e+14 -8.41037e+14 -8.42169e+14 -8.43298e+14 -8.44425e+14 -8.4555e+14 -8.46672e+14 -8.47792e+14 -8.48911e+14 -8.50026e+14 -8.5114e+14 -8.52251e+14 -8.5336e+14 -8.54467e+14 -8.55571e+14 -8.56674e+14 -8.57774e+14 -8.58871e+14 -8.59967e+14 -8.6106e+14 -8.62151e+14 -8.6324e+14 -8.64326e+14 -8.65411e+14 -8.66493e+14 -8.67572e+14 -8.6865e+14 -8.69725e+14 -8.70798e+14 -8.71869e+14 -8.72937e+14 -8.74003e+14 -8.75067e+14 -8.76129e+14 -8.77189e+14 -8.78246e+14 -8.79301e+14 -8.80354e+14 -8.81404e+14 -8.82452e+14 -8.83498e+14 -8.84542e+14 -8.85583e+14 -8.86623e+14 -8.87659e+14 -8.88694e+14 -8.89727e+14 -8.90757e+14 -8.91785e+14 -8.9281e+14 -8.93834e+14 -8.94855e+14 -8.95874e+14 -8.96891e+14 -8.97905e+14 -8.98917e+14 -8.99927e+14 -9.00935e+14 -9.0194e+14 -9.02944e+14 -9.03945e+14 -9.04943e+14 -9.0594e+14 -9.06934e+14 -9.07926e+14 -9.08915e+14 -9.09903e+14 -9.10888e+14 -9.11871e+14 -9.12851e+14 -9.1383e+14 -9.14806e+14 -9.1578e+14 -9.16752e+14 -9.17721e+14 -9.18688e+14 -9.19653e+14 -9.20616e+14 -9.21576e+14 -9.22534e+14 -9.2349e+14 -9.24444e+14 -9.25395e+14 -9.26344e+14 -9.27291e+14 -9.28236e+14 -9.29178e+14 -9.30118e+14 -9.31056e+14 -9.31992e+14 -9.32925e+14 -9.33856e+14 -9.34785e+14 -9.35712e+14 -9.36636e+14 -9.37558e+14 -9.38478e+14 -9.39396e+14 -9.40311e+14 -9.41224e+14 -9.42135e+14 -9.43044e+14 -9.4395e+14 -9.44854e+14 -9.45756e+14 -9.46656e+14 -9.47553e+14 -9.48448e+14 -9.49341e+14 -9.50231e+14 -9.5112e+14 -9.52006e+14 -9.5289e+14 -9.53771e+14 -9.54651e+14 -9.55528e+14 -9.56402e+14 -9.57275e+14 -9.58145e+14 -9.59013e+14 -9.59879e+14 -9.60743e+14 -9.61604e+14 -9.62463e+14 -9.6332e+14 -9.64175e+14 -9.65027e+14 -9.65877e+14 -9.66725e+14 -9.6757e+14 -9.68414e+14 -9.69255e+14 -9.70093e+14 -9.7093e+14 -9.71764e+14 -9.72596e+14 -9.73426e+14 -9.74254e+14 -9.75079e+14 -9.75902e+14 -9.76723e+14 -9.77541e+14 -9.78357e+14 -9.79172e+14 -9.79983e+14 -9.80793e+14 -9.816e+14 -9.82405e+14 -9.83208e+14 -9.84008e+14 -9.84807e+14 -9.85603e+14 -9.86396e+14 -9.87188e+14 -9.87977e+14 -9.88764e+14 -9.89549e+14 -9.90331e+14 -9.91112e+14 -9.9189e+14 -9.92665e+14 -9.93439e+14 -9.9421e+14 -9.94979e+14 -9.95746e+14 -9.9651e+14 -9.97273e+14 -9.98033e+14 -9.9879e+14 -9.99546e+14 -1.0003e+15 -1.00105e+15 -1.0018e+15 -1.00255e+15 -1.00329e+15 -1.00403e+15 -1.00477e+15 -1.00551e+15 -1.00624e+15 -1.00698e+15 -1.00771e+15 -1.00844e+15 -1.00916e+15 -1.00989e+15 -1.01061e+15 -1.01133e+15 -1.01204e+15 -1.01276e+15 -1.01347e+15 -1.01418e+15 -1.01489e+15 -1.0156e+15 -1.0163e+15 -1.017e+15 -1.0177e+15 -1.0184e+15 -1.01909e+15 -1.01978e+15 -1.02048e+15 -1.02116e+15 -1.02185e+15 -1.02253e+15 -1.02321e+15 -1.02389e+15 -1.02457e+15 -1.02524e+15 -1.02592e+15 -1.02659e+15 -1.02725e+15 -1.02792e+15 -1.02858e+15 -1.02924e+15 -1.0299e+15 -1.03056e+15 -1.03121e+15 -1.03186e+15 -1.03251e+15 -1.03316e+15 -1.03381e+15 -1.03445e+15 -1.03509e+15 -1.03573e+15 -1.03636e+15 -1.037e+15 -1.03763e+15 -1.03826e+15 -1.03889e+15 -1.03951e+15 -1.04013e+15 -1.04075e+15 -1.04137e+15 -1.04199e+15 -1.0426e+15 -1.04321e+15 -1.04382e+15 -1.04443e+15 -1.04503e+15 -1.04564e+15 -1.04624e+15 -1.04683e+15 -1.04743e+15 -1.04802e+15 -1.04861e+15 -1.0492e+15 -1.04979e+15 -1.05037e+15 -1.05096e+15 -1.05154e+15 -1.05211e+15 -1.05269e+15 -1.05326e+15 -1.05383e+15 -1.0544e+15 -1.05497e+15 -1.05553e+15 -1.05609e+15 -1.05665e+15 -1.05721e+15 -1.05776e+15 -1.05832e+15 -1.05887e+15 -1.05942e+15 -1.05996e+15 -1.06051e+15 -1.06105e+15 -1.06159e+15 -1.06212e+15 -1.06266e+15 -1.06319e+15 -1.06372e+15 -1.06425e+15 -1.06478e+15 -1.0653e+15 -1.06582e+15 -1.06634e+15 -1.06686e+15 -1.06737e+15 -1.06788e+15 -1.06839e+15 -1.0689e+15 -1.06941e+15 -1.06991e+15 -1.07041e+15 -1.07091e+15 -1.07141e+15 -1.0719e+15 -1.07239e+15 -1.07288e+15 -1.07337e+15 -1.07386e+15 -1.07434e+15 -1.07482e+15 -1.0753e+15 -1.07577e+15 -1.07625e+15 -1.07672e+15 -1.07719e+15 -1.07766e+15 -1.07812e+15 -1.07858e+15 -1.07904e+15 -1.0795e+15 -1.07996e+15 -1.08041e+15 -1.08086e+15 -1.08131e+15 -1.08176e+15 -1.0822e+15 -1.08265e+15 -1.08309e+15 -1.08353e+15 -1.08396e+15 -1.08439e+15 -1.08483e+15 -1.08525e+15 -1.08568e+15 -1.08611e+15 -1.08653e+15 -1.08695e+15 -1.08737e+15 -1.08778e+15 -1.08819e+15 -1.08861e+15 -1.08901e+15 -1.08942e+15 -1.08982e+15 -1.09023e+15 -1.09063e+15 -1.09102e+15 -1.09142e+15 -1.09181e+15 -1.0922e+15 -1.09259e+15 -1.09298e+15 -1.09336e+15 -1.09374e+15 -1.09412e+15 -1.0945e+15 -1.09487e+15 -1.09525e+15 -1.09562e+15 -1.09599e+15 -1.09635e+15 -1.09672e+15 -1.09708e+15 -1.09744e+15 -1.09779e+15 -1.09815e+15 -1.0985e+15 -1.09885e+15 -1.0992e+15 -1.09954e+15 -1.09989e+15 -1.10023e+15 -1.10057e+15 -1.1009e+15 -1.10124e+15 -1.10157e+15 -1.1019e+15 -1.10223e+15 -1.10255e+15 -1.10288e+15 -1.1032e+15 -1.10352e+15 -1.10383e+15 -1.10415e+15 -1.10446e+15 -1.10477e+15 -1.10508e+15 -1.10538e+15 -1.10568e+15 -1.10598e+15 -1.10628e+15 -1.10658e+15 -1.10687e+15 -1.10716e+15 -1.10745e+15 -1.10774e+15 -1.10803e+15 -1.10831e+15 -1.10859e+15 -1.10887e+15 -1.10914e+15 -1.10942e+15 -1.10969e+15 -1.10996e+15 -1.11022e+15 -1.11049e+15 -1.11075e+15 -1.11101e+15 -1.11127e+15 -1.11152e+15 -1.11178e+15 -1.11203e+15 -1.11228e+15 -1.11252e+15 -1.11277e+15 -1.11301e+15 -1.11325e+15 -1.11349e+15 -1.11372e+15 -1.11396e+15 -1.11419e+15 -1.11441e+15 -1.11464e+15 -1.11486e+15 -1.11509e+15 -1.11531e+15 -1.11552e+15 -1.11574e+15 -1.11595e+15 -1.11616e+15 -1.11637e+15 -1.11658e+15 -1.11678e+15 -1.11698e+15 -1.11718e+15 -1.11738e+15 -1.11757e+15 -1.11777e+15 -1.11796e+15 -1.11814e+15 -1.11833e+15 -1.11851e+15 -1.11869e+15 -1.11887e+15 -1.11905e+15 -1.11922e+15 -1.1194e+15 -1.11957e+15 -1.11973e+15 -1.1199e+15 -1.12006e+15 -1.12022e+15 -1.12038e+15 -1.12054e+15 -1.12069e+15 -1.12085e+15 -1.121e+15 -1.12114e+15 -1.12129e+15 -1.12143e+15 -1.12157e+15 -1.12171e+15 -1.12185e+15 -1.12198e+15 -1.12211e+15 -1.12224e+15 -1.12237e+15 -1.12249e+15 -1.12262e+15 -1.12274e+15 -1.12286e+15 -1.12297e+15 -1.12309e+15 -1.1232e+15 -1.12331e+15 -1.12341e+15 -1.12352e+15 -1.12362e+15 -1.12372e+15 -1.12382e+15 -1.12391e+15 -1.12401e+15 -1.1241e+15 -1.12419e+15 -1.12427e+15 -1.12436e+15 -1.12444e+15 -1.12452e+15 -1.1246e+15 -1.12467e+15 -1.12475e+15 -1.12482e+15 -1.12489e+15 -1.12495e+15 -1.12502e+15 -1.12508e+15 -1.12514e+15 -1.1252e+15 -1.12525e+15 -1.1253e+15 -1.12535e+15 -1.1254e+15 -1.12545e+15 -1.12549e+15 -1.12554e+15 -1.12557e+15 -1.12561e+15 -1.12565e+15 -1.12568e+15 -1.12571e+15 -1.12574e+15 -1.12576e+15 -1.12579e+15 -1.12581e+15 -1.12583e+15 -1.12584e+15 -1.12586e+15 -1.12587e+15 -1.12588e+15 -1.12589e+15 -1.1259e+15 -1.1259e+15 -1.1259e+15 -1.1259e+15 -1.1259e+15 -1.12589e+15 -1.12588e+15 -1.12587e+15 -1.12586e+15 -1.12584e+15 -1.12583e+15 -1.12581e+15 -1.12579e+15 -1.12576e+15 -1.12574e+15 -1.12571e+15 -1.12568e+15 -1.12565e+15 -1.12561e+15 -1.12557e+15 -1.12554e+15 -1.12549e+15 -1.12545e+15 -1.1254e+15 -1.12535e+15 -1.1253e+15 -1.12525e+15 -1.1252e+15 -1.12514e+15 -1.12508e+15 -1.12502e+15 -1.12495e+15 -1.12489e+15 -1.12482e+15 -1.12475e+15 -1.12467e+15 -1.1246e+15 -1.12452e+15 -1.12444e+15 -1.12436e+15 -1.12427e+15 -1.12419e+15 -1.1241e+15 -1.12401e+15 -1.12391e+15 -1.12382e+15 -1.12372e+15 -1.12362e+15 -1.12352e+15 -1.12341e+15 -1.12331e+15 -1.1232e+15 -1.12309e+15 -1.12297e+15 -1.12286e+15 -1.12274e+15 -1.12262e+15 -1.12249e+15 -1.12237e+15 -1.12224e+15 -1.12211e+15 -1.12198e+15 -1.12185e+15 -1.12171e+15 -1.12157e+15 -1.12143e+15 -1.12129e+15 -1.12114e+15 -1.121e+15 -1.12085e+15 -1.12069e+15 -1.12054e+15 -1.12038e+15 -1.12022e+15 -1.12006e+15 -1.1199e+15 -1.11973e+15 -1.11957e+15 -1.1194e+15 -1.11922e+15 -1.11905e+15 -1.11887e+15 -1.11869e+15 -1.11851e+15 -1.11833e+15 -1.11814e+15 -1.11796e+15 -1.11777e+15 -1.11757e+15 -1.11738e+15 -1.11718e+15 -1.11698e+15 -1.11678e+15 -1.11658e+15 -1.11637e+15 -1.11616e+15 -1.11595e+15 -1.11574e+15 -1.11552e+15 -1.11531e+15 -1.11509e+15 -1.11486e+15 -1.11464e+15 -1.11441e+15 -1.11419e+15 -1.11396e+15 -1.11372e+15 -1.11349e+15 -1.11325e+15 -1.11301e+15 -1.11277e+15 -1.11252e+15 -1.11228e+15 -1.11203e+15 -1.11178e+15 -1.11152e+15 -1.11127e+15 -1.11101e+15 -1.11075e+15 -1.11049e+15 -1.11022e+15 -1.10996e+15 -1.10969e+15 -1.10942e+15 -1.10914e+15 -1.10887e+15 -1.10859e+15 -1.10831e+15 -1.10803e+15 -1.10774e+15 -1.10745e+15 -1.10716e+15 -1.10687e+15 -1.10658e+15 -1.10628e+15 -1.10598e+15 -1.10568e+15 -1.10538e+15 -1.10508e+15 -1.10477e+15 -1.10446e+15 -1.10415e+15 -1.10383e+15 -1.10352e+15 -1.1032e+15 -1.10288e+15 -1.10255e+15 -1.10223e+15 -1.1019e+15 -1.10157e+15 -1.10124e+15 -1.1009e+15 -1.10057e+15 -1.10023e+15 -1.09989e+15 -1.09954e+15 -1.0992e+15 -1.09885e+15 -1.0985e+15 -1.09815e+15 -1.09779e+15 -1.09744e+15 -1.09708e+15 -1.09672e+15 -1.09635e+15 -1.09599e+15 -1.09562e+15 -1.09525e+15 -1.09487e+15 -1.0945e+15 -1.09412e+15 -1.09374e+15 -1.09336e+15 -1.09298e+15 -1.09259e+15 -1.0922e+15 -1.09181e+15 -1.09142e+15 -1.09102e+15 -1.09063e+15 -1.09023e+15 -1.08982e+15 -1.08942e+15 -1.08901e+15 -1.08861e+15 -1.08819e+15 -1.08778e+15 -1.08737e+15 -1.08695e+15 -1.08653e+15 -1.08611e+15 -1.08568e+15 -1.08525e+15 -1.08483e+15 -1.08439e+15 -1.08396e+15 -1.08353e+15 -1.08309e+15 -1.08265e+15 -1.0822e+15 -1.08176e+15 -1.08131e+15 -1.08086e+15 -1.08041e+15 -1.07996e+15 -1.0795e+15 -1.07904e+15 -1.07858e+15 -1.07812e+15 -1.07766e+15 -1.07719e+15 -1.07672e+15 -1.07625e+15 -1.07577e+15 -1.0753e+15 -1.07482e+15 -1.07434e+15 -1.07386e+15 -1.07337e+15 -1.07288e+15 -1.07239e+15 -1.0719e+15 -1.07141e+15 -1.07091e+15 -1.07041e+15 -1.06991e+15 -1.06941e+15 -1.0689e+15 -1.06839e+15 -1.06788e+15 -1.06737e+15 -1.06686e+15 -1.06634e+15 -1.06582e+15 -1.0653e+15 -1.06478e+15 -1.06425e+15 -1.06372e+15 -1.06319e+15 -1.06266e+15 -1.06212e+15 -1.06159e+15 -1.06105e+15 -1.06051e+15 -1.05996e+15 -1.05942e+15 -1.05887e+15 -1.05832e+15 -1.05776e+15 -1.05721e+15 -1.05665e+15 -1.05609e+15 -1.05553e+15 -1.05497e+15 -1.0544e+15 -1.05383e+15 -1.05326e+15 -1.05269e+15 -1.05211e+15 -1.05154e+15 -1.05096e+15 -1.05037e+15 -1.04979e+15 -1.0492e+15 -1.04861e+15 -1.04802e+15 -1.04743e+15 -1.04683e+15 -1.04624e+15 -1.04564e+15 -1.04503e+15 -1.04443e+15 -1.04382e+15 -1.04321e+15 -1.0426e+15 -1.04199e+15 -1.04137e+15 -1.04075e+15 -1.04013e+15 -1.03951e+15 -1.03889e+15 -1.03826e+15 -1.03763e+15 -1.037e+15 -1.03636e+15 -1.03573e+15 -1.03509e+15 -1.03445e+15 -1.03381e+15 -1.03316e+15 -1.03251e+15 -1.03186e+15 -1.03121e+15 -1.03056e+15 -1.0299e+15 -1.02924e+15 -1.02858e+15 -1.02792e+15 -1.02725e+15 -1.02659e+15 -1.02592e+15 -1.02524e+15 -1.02457e+15 -1.02389e+15 -1.02321e+15 -1.02253e+15 -1.02185e+15 -1.02116e+15 -1.02048e+15 -1.01978e+15 -1.01909e+15 -1.0184e+15 -1.0177e+15 -1.017e+15 -1.0163e+15 -1.0156e+15 -1.01489e+15 -1.01418e+15 -1.01347e+15 -1.01276e+15 -1.01204e+15 -1.01133e+15 -1.01061e+15 -1.00989e+15 -1.00916e+15 -1.00844e+15 -1.00771e+15 -1.00698e+15 -1.00624e+15 -1.00551e+15 -1.00477e+15 -1.00403e+15 -1.00329e+15 -1.00255e+15 -1.0018e+15 -1.00105e+15 -1.0003e+15 -9.99546e+14 -9.9879e+14 -9.98033e+14 -9.97273e+14 -9.9651e+14 -9.95746e+14 -9.94979e+14 -9.9421e+14 -9.93439e+14 -9.92665e+14 -9.9189e+14 -9.91112e+14 -9.90331e+14 -9.89549e+14 -9.88764e+14 -9.87977e+14 -9.87188e+14 -9.86396e+14 -9.85603e+14 -9.84807e+14 -9.84008e+14 -9.83208e+14 -9.82405e+14 -9.816e+14 -9.80793e+14 -9.79983e+14 -9.79172e+14 -9.78357e+14 -9.77541e+14 -9.76723e+14 -9.75902e+14 -9.75079e+14 -9.74254e+14 -9.73426e+14 -9.72596e+14 -9.71764e+14 -9.7093e+14 -9.70093e+14 -9.69255e+14 -9.68414e+14 -9.6757e+14 -9.66725e+14 -9.65877e+14 -9.65027e+14 -9.64175e+14 -9.6332e+14 -9.62463e+14 -9.61604e+14 -9.60743e+14 -9.59879e+14 -9.59013e+14 -9.58145e+14 -9.57275e+14 -9.56402e+14 -9.55528e+14 -9.54651e+14 -9.53771e+14 -9.5289e+14 -9.52006e+14 -9.5112e+14 -9.50231e+14 -9.49341e+14 -9.48448e+14 -9.47553e+14 -9.46656e+14 -9.45756e+14 -9.44854e+14 -9.4395e+14 -9.43044e+14 -9.42135e+14 -9.41224e+14 -9.40311e+14 -9.39396e+14 -9.38478e+14 -9.37558e+14 -9.36636e+14 -9.35712e+14 -9.34785e+14 -9.33856e+14 -9.32925e+14 -9.31992e+14 -9.31056e+14 -9.30118e+14 -9.29178e+14 -9.28236e+14 -9.27291e+14 -9.26344e+14 -9.25395e+14 -9.24444e+14 -9.2349e+14 -9.22534e+14 -9.21576e+14 -9.20616e+14 -9.19653e+14 -9.18688e+14 -9.17721e+14 -9.16752e+14 -9.1578e+14 -9.14806e+14 -9.1383e+14 -9.12851e+14 -9.11871e+14 -9.10888e+14 -9.09903e+14 -9.08915e+14 -9.07926e+14 -9.06934e+14 -9.0594e+14 -9.04943e+14 -9.03945e+14 -9.02944e+14 -9.0194e+14 -9.00935e+14 -8.99927e+14 -8.98917e+14 -8.97905e+14 -8.96891e+14 -8.95874e+14 -8.94855e+14 -8.93834e+14 -8.9281e+14 -8.91785e+14 -8.90757e+14 -8.89727e+14 -8.88694e+14 -8.87659e+14 -8.86623e+14 -8.85583e+14 -8.84542e+14 -8.83498e+14 -8.82452e+14 -8.81404e+14 -8.80354e+14 -8.79301e+14 -8.78246e+14 -8.77189e+14 -8.76129e+14 -8.75067e+14 -8.74003e+14 -8.72937e+14 -8.71869e+14 -8.70798e+14 -8.69725e+14 -8.6865e+14 -8.67572e+14 -8.66493e+14 -8.65411e+14 -8.64326e+14 -8.6324e+14 -8.62151e+14 -8.6106e+14 -8.59967e+14 -8.58871e+14 -8.57774e+14 -8.56674e+14 -8.55571e+14 -8.54467e+14 -8.5336e+14 -8.52251e+14 -8.5114e+14 -8.50026e+14 -8.48911e+14 -8.47792e+14 -8.46672e+14 -8.4555e+14 -8.44425e+14 -8.43298e+14 -8.42169e+14 -8.41037e+14 -8.39903e+14 -8.38767e+14 -8.37629e+14 -8.36488e+14 -8.35346e+14 -8.34201e+14 -8.33053e+14 -8.31904e+14 -8.30752e+14 -8.29598e+14 -8.28442e+14 -8.27283e+14 -8.26122e+14 -8.24959e+14 -8.23794e+14 -8.22626e+14 -8.21457e+14 -8.20285e+14 -8.1911e+14 -8.17934e+14 -8.16755e+14 -8.15574e+14 -8.1439e+14 -8.13205e+14 -8.12017e+14 -8.10827e+14 -8.09635e+14 -8.0844e+14 -8.07243e+14 -8.06044e+14 -8.04843e+14 -8.03639e+14 -8.02433e+14 -8.01225e+14 -8.00015e+14 -7.98802e+14 -7.97587e+14 -7.9637e+14 -7.95151e+14 -7.93929e+14 -7.92706e+14 -7.91479e+14 -7.90251e+14 -7.89021e+14 -7.87788e+14 -7.86553e+14 -7.85315e+14 -7.84076e+14 -7.82834e+14 -7.8159e+14 -7.80343e+14 -7.79095e+14 -7.77844e+14 -7.76591e+14 -7.75335e+14 -7.74078e+14 -7.72818e+14 -7.71556e+14 -7.70291e+14 -7.69025e+14 -7.67756e+14 -7.66485e+14 -7.65211e+14 -7.63935e+14 -7.62658e+14 -7.61377e+14 -7.60095e+14 -7.5881e+14 -7.57523e+14 -7.56234e+14 -7.54943e+14 -7.53649e+14 -7.52353e+14 -7.51055e+14 -7.49755e+14 -7.48452e+14 -7.47147e+14 -7.4584e+14 -7.44531e+14 -7.43219e+14 -7.41905e+14 -7.40589e+14 -7.3927e+14 -7.3795e+14 -7.36627e+14 -7.35302e+14 -7.33974e+14 -7.32644e+14 -7.31313e+14 -7.29978e+14 -7.28642e+14 -7.27303e+14 -7.25962e+14 -7.24619e+14 -7.23274e+14 -7.21926e+14 -7.20576e+14 -7.19224e+14 -7.17869e+14 -7.16513e+14 -7.15154e+14 -7.13792e+14 -7.12429e+14 -7.11063e+14 -7.09695e+14 -7.08325e+14 -7.06953e+14 -7.05578e+14 -7.04201e+14 -7.02822e+14 -7.0144e+14 -7.00056e+14 -6.9867e+14 -6.97282e+14 -6.95892e+14 -6.94499e+14 -6.93104e+14 -6.91707e+14 -6.90307e+14 -6.88906e+14 -6.87502e+14 -6.86095e+14 -6.84687e+14 -6.83276e+14 -6.81863e+14 -6.80448e+14 -6.7903e+14 -6.7761e+14 -6.76188e+14 -6.74764e+14 -6.73338e+14 -6.71909e+14 -6.70478e+14 -6.69045e+14 -6.67609e+14 -6.66171e+14 -6.64731e+14 -6.63289e+14 -6.61844e+14 -6.60398e+14 -6.58949e+14 -6.57497e+14 -6.56044e+14 -6.54588e+14 -6.5313e+14 -6.5167e+14 -6.50207e+14 -6.48742e+14 -6.47275e+14 -6.45806e+14 -6.44335e+14 -6.42861e+14 -6.41385e+14 -6.39906e+14 -6.38426e+14 -6.36943e+14 -6.35458e+14 -6.33971e+14 -6.32481e+14 -6.30989e+14 -6.29495e+14 -6.27999e+14 -6.265e+14 -6.24999e+14 -6.23496e+14 -6.21991e+14 -6.20483e+14 -6.18974e+14 -6.17462e+14 -6.15947e+14 -6.14431e+14 -6.12912e+14 -6.11391e+14 -6.09867e+14 -6.08342e+14 -6.06814e+14 -6.05284e+14 -6.03751e+14 -6.02217e+14 -6.0068e+14 -5.99141e+14 -5.976e+14 -5.96056e+14 -5.9451e+14 -5.92962e+14 -5.91412e+14 -5.89859e+14 -5.88304e+14 -5.86747e+14 -5.85188e+14 -5.83626e+14 -5.82062e+14 -5.80496e+14 -5.78928e+14 -5.77357e+14 -5.75784e+14 -5.74209e+14 -5.72632e+14 -5.71052e+14 -5.6947e+14 -5.67886e+14 -5.663e+14 -5.64711e+14 -5.6312e+14 -5.61527e+14 -5.59931e+14 -5.58334e+14 -5.56734e+14 -5.55132e+14 -5.53527e+14 -5.51921e+14 -5.50312e+14 -5.48701e+14 -5.47087e+14 -5.45471e+14 -5.43854e+14 -5.42233e+14 -5.40611e+14 -5.38986e+14 -5.37359e+14 -5.3573e+14 -5.34099e+14 -5.32465e+14 -5.30829e+14 -5.29191e+14 -5.27551e+14 -5.25908e+14 -5.24263e+14 -5.22616e+14 -5.20966e+14 -5.19315e+14 -5.17661e+14 -5.16004e+14 -5.14346e+14 -5.12685e+14 -5.11022e+14 -5.09357e+14 -5.0769e+14 -5.0602e+14 -5.04348e+14 -5.02674e+14 -5.00997e+14 -4.99319e+14 -4.97638e+14 -4.95954e+14 -4.94269e+14 -4.92581e+14 -4.90891e+14 -4.89199e+14 -4.87505e+14 -4.85808e+14 -4.84109e+14 -4.82408e+14 -4.80704e+14 -4.78998e+14 -4.7729e+14 -4.7558e+14 -4.73868e+14 -4.72153e+14 -4.70436e+14 -4.68717e+14 -4.66995e+14 -4.65271e+14 -4.63545e+14 -4.61817e+14 -4.60087e+14 -4.58354e+14 -4.56619e+14 -4.54882e+14 -4.53142e+14 -4.514e+14 -4.49656e+14 -4.4791e+14 -4.46161e+14 -4.44411e+14 -4.42658e+14 -4.40902e+14 -4.39145e+14 -4.37385e+14 -4.35623e+14 -4.33859e+14 -4.32092e+14 -4.30323e+14 -4.28552e+14 -4.26779e+14 -4.25004e+14 -4.23226e+14 -4.21446e+14 -4.19663e+14 -4.17879e+14 -4.16092e+14 -4.14303e+14 -4.12512e+14 -4.10718e+14 -4.08922e+14 -4.07124e+14 -4.05324e+14 -4.03521e+14 -4.01717e+14 -3.9991e+14 -3.981e+14 -3.96289e+14 -3.94475e+14 -3.92659e+14 -3.9084e+14 -3.8902e+14 -3.87197e+14 -3.85372e+14 -3.83545e+14 -3.81715e+14 -3.79883e+14 -3.78049e+14 -3.76213e+14 -3.74374e+14 -3.72533e+14 -3.7069e+14 -3.68845e+14 -3.66997e+14 -3.65147e+14 -3.63295e+14 -3.61441e+14 -3.59584e+14 -3.57725e+14 -3.55864e+14 -3.54001e+14 -3.52135e+14 -3.50267e+14 -3.48397e+14 -3.46525e+14 -3.4465e+14 -3.42773e+14 -3.40894e+14 -3.39013e+14 -3.37129e+14 -3.35243e+14 -3.33355e+14 -3.31465e+14 -3.29572e+14 -3.27677e+14 -3.2578e+14 -3.23881e+14 -3.21979e+14 -3.20075e+14 -3.18169e+14 -3.16261e+14 -3.1435e+14 -3.12437e+14 -3.10522e+14 -3.08605e+14 -3.06685e+14 -3.04763e+14 -3.02839e+14 -3.00913e+14 -2.98984e+14 -2.97053e+14 -2.9512e+14 -2.93184e+14 -2.91247e+14 -2.89307e+14 -2.87365e+14 -2.8542e+14 -2.83473e+14 -2.81525e+14 -2.79573e+14 -2.7762e+14 -2.75664e+14 -2.73706e+14 -2.71746e+14 -2.69784e+14 -2.67819e+14 -2.65852e+14 -2.63883e+14 -2.61911e+14 -2.59938e+14 -2.57962e+14 -2.55983e+14 -2.54003e+14 -2.5202e+14 -2.50035e+14 -2.48048e+14 -2.46059e+14 -2.44067e+14 -2.42073e+14 -2.40077e+14 -2.38078e+14 -2.36078e+14 -2.34075e+14 -2.32069e+14 -2.30062e+14 -2.28052e+14 -2.2604e+14 -2.24026e+14 -2.22009e+14 -2.19991e+14 -2.1797e+14 -2.15946e+14 -2.13921e+14 -2.11893e+14 -2.09863e+14 -2.07831e+14 -2.05796e+14 -2.0376e+14 -2.01721e+14 -1.99679e+14 -1.97636e+14 -1.9559e+14 -1.93542e+14 -1.91492e+14 -1.89439e+14 -1.87385e+14 -1.85328e+14 -1.83268e+14 -1.81207e+14 -1.79143e+14 -1.77077e+14 -1.75009e+14 -1.72938e+14 -1.70865e+14 -1.6879e+14 -1.66713e+14 -1.64634e+14 -1.62552e+14 -1.60468e+14 -1.58381e+14 -1.56293e+14 -1.54202e+14 -1.52109e+14 -1.50014e+14 -1.47916e+14 -1.45816e+14 -1.43714e+14 -1.4161e+14 -1.39504e+14 -1.37395e+14 -1.35284e+14 -1.3317e+14 -1.31055e+14 -1.28937e+14 -1.26817e+14 -1.24695e+14 -1.2257e+14 -1.20443e+14 -1.18314e+14 -1.16183e+14 -1.14049e+14 -1.11913e+14 -1.09775e+14 -1.07635e+14 -1.05492e+14 -1.03347e+14 -1.012e+14 -9.9051e+13 -9.68994e+13 -9.47456e+13 -9.25895e+13 -9.04312e+13 -8.82706e+13 -8.61077e+13 -8.39426e+13 -8.17752e+13 -7.96056e+13 -7.74338e+13 -7.52597e+13 -7.30833e+13 -7.09047e+13 -6.87238e+13 -6.65407e+13 -6.43553e+13 -6.21677e+13 -5.99778e+13 -5.77857e+13 -5.55913e+13 -5.33947e+13 -5.11958e+13 -4.89947e+13 -4.67913e+13 -4.45856e+13 -4.23777e+13 -4.01676e+13 -3.79552e+13 -3.57406e+13 -3.35237e+13 -3.13045e+13 -2.90831e+13 -2.68595e+13 -2.46336e+13 -2.24054e+13 -2.0175e+13 -1.79423e+13 -1.57074e+13 -1.34703e+13 -1.12309e+13 -8.98918e+12 -6.74527e+12 -4.4991e+12 -2.25067e+12]
plot(x,z);

%%

clear all;
close all;
clc;

uex=@(x,y)[sin(2*pi*x).*sin(4*pi*y)]; % exact solution = boundary data 
uex_x=@(x,y)[2 * pi *cos(2*pi*x).*sin(4*pi*y)]; % du/dx 
uex_y=@(x,y)[4 * pi * sin(2*pi*x).*cos(4*pi*y)]; % du/dy 
f=@(x,y)[(20*pi*pi + 1)*sin(2*pi*x).*sin(4*pi*y)]; % r.h.s 
g=@(x,y)[0]; % Dirichlet boundary data 
h=@(x,y)[0]; % Neumann boundary data 
gam=1; % coefficient of the term of order zero.
xa=0;xb=1; % Omega=(xa,xb) x (ya,yb)
ya=0;yb=1;
cb='dddd'; % For boundary conditions cb(i)=;d’ ---> Dirichlet on side i
% cb(i)=’n’ ----> Neumann on side i 
nex=10;ney=10; nx=5; ny=5;
param=zeros(20,1); % "help lap_2d" for a complete description of param
param(1)=1;
param(2)=2;
param(3)=1;
% 1=SEM-NI, 2= Patching
% 0=no reordering, 1=CM ordering, 2=AMD ordering
% 1= solve linear system by Cholesky fact.
% 2= compute extrema eigenvalues of A
% 3 solve by Schur complement
% 4= compute extrema eigenvalues of the Schur complement
param(4)=1; % computes errors
param(5)=1; % 0 exact norms, 1= discrete norms
param(6)=nx*2;
param(7)=1; %
param(8)=2; %
param(9)=(nx+1); % nodes used to plot numerical solution 
gammax=[]; gammay=[]; % if SEM decomposition is not uniform:
% nq for LG quadrature formulas
%0 =absolute errors, 1=relative errors 0 no plot, 1 mesh, 2 surf, 3 contour
% they are the arrays with intefaces positions % along x- and y- directions, respectively


[xy,un,D,param, A, F, weights]=lap_2d(xa,xb,ya,yb,gam,uex,uex_x,uex_y,f,g,h,cb,nex,nx,ney,ny,gammax,gammay,param);

fprintf('nx=%d,nex=%d,err_inf=%11.4e, err_h1=%11.4e,err_l2=%11.4e \n', nx,nex,param(29),param(30),param(31))

A = full(A);

%write the matrix to a file to compare with c++ matrix

filename = 'A.csv';
csvwrite(filename, A);
filename = 'F.csv';
csvwrite(filename, F);

%%

 x = [1,2,4,7];
 x = reshape(x,1,1)
