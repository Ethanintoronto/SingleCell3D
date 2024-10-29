clear all;
%This script goes through the derivation of the forces in a 3D vertex
%model.

% We start with the energy functional E = K_V(V-V_0)^2 + K_A(A-A_0)^2 
% We will then solve:
% F = -grad(E) --> F = - (2K_V (V-V_0)*grad(V) + 2K_A(A-A_0)*grad(A))
% Therefore we now write V and A as functions of the cartesian coordinates
syms x_i y_i z_i x_j y_j z_j Cx_0 Cy_0 Cz_0 Px_0 Py_0 Pz_0 N_p x_k y_k z_k real
r_i = [x_i,y_i,z_i];
r_j = [x_j,y_j, z_j];
C_0 = [Cx_0, Cy_0, Cz_0];
P_0 = [Px_0, Py_0, Pz_0]; 
%let [x_i,y_i,z_i] be a general vertex
%let [x_j,y_j,z_j] be the next counterclockwise vertex after i such that
%j = i+1%nV where nV is the total number of vertices on the face

%let Cx_0 Cy_0 Cz_0 be the reference coordinates for cell:
%The first vertex of the first polygon - note that we do not require a
%del_Ck since the volume contribution is zero if vertex i or j is the
%reference vertex

%let Px_0 Py_0 Pz_0 be the reference coordinates for the polygon:
%The centroid of given face

%Let N_p be the number of vertices on a given face p.

% VOLUME DERIVATION:
% We can write the total cell volume as 
% the sum of the absolute value of the triple product of all tetrahedrons
% formed by the vector from reference points V_0 to A_0 
% and the vertices of each edge ([x_i,y_i,z_i] , [x_j,y_j,z_j])
V_tetra = dot(cross((r_i-P_0), (r_j-P_0)), (P_0-C_0));

% Collecting the volume equation in terms of x_i and x_j:
V_tetra_x = collect(V_tetra,[x_i,x_j,Px_0]);
disp(V_tetra_x)

V_terms_x = children(V_tetra_x); %term 1 = x_i, term 2 = x_j, term 3 = Px_0 term 4 = none

% Apply the derivative 
dVdx_k = kronnecker_delta(diff(V_terms_x{1},x_i), {'k', 'kp1'});
dVdx_k = dVdx_k + kronnecker_delta(diff(V_terms_x{2},x_j), {'km1','k'});
dVdx_k = dVdx_k + kronnecker_delta(diff(V_terms_x{3},Px_0), {'k','kp1'})/N_p;
% AREA DERIVATION:
% We can write the area of the cell as the sum of the areas of each polygon
% The faces of each polygon can be written as the sum of the areas of the 
% triangles formed by the centroid and all edge vertices.

A_tri = sum(cross((r_i-P_0), (r_j-P_0)).^2); 
A_tri_x = collect(expand(A_tri), [x_i,x_j, x_i^2, x_j^2,Px_0, Px_0^2]);

A_terms_x = children(A_tri_x); %term 1 = x_i*x_j term2 = x_i*Px_0 term3 = x_jPx_0 term4 = x_i^2 term5 = x_j^2 term 6 = Px_0^2 term 7 = no x

dAdx_k = kronnecker_delta(diff(A_terms_x{1},x_i),{'k', 'kp1'}); 
dAdx_k = dAdx_k + kronnecker_delta(diff(A_terms_x{1},x_j),{'km1', 'k'});

dAdx_k = dAdx_k + kronnecker_delta(diff(A_terms_x{2},x_i), {'k', 'kp1'});
dAdx_k = dAdx_k + kronnecker_delta(diff(A_terms_x{2},Px_0), {'k', 'kp1'})/N_p;

dAdx_k = dAdx_k + kronnecker_delta(diff(A_terms_x{3},x_j), {'km1', 'k'});
dAdx_k = dAdx_k + kronnecker_delta(diff(A_terms_x{3},Px_0), {'k', 'kp1'})/N_p;

dAdx_k = dAdx_k + kronnecker_delta(diff(A_terms_x{4},x_i), {'k', 'kp1'});

dAdx_k = dAdx_k + kronnecker_delta(diff(A_terms_x{5},x_j), {'km1', 'k'});

dAdx_k = dAdx_k + kronnecker_delta(diff(A_terms_x{6},Px_0), {'k', 'kp1'})/N_p;

dAdx_k = 1/4*(dAdx_k)*A_tri_x^(-1/2); %the 1/4 comes from 1/2 the area to get the triangle and 1/2 from taking the derivative



