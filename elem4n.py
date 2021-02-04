# -*- coding: utf-8 -*-
"""
Created on Thu Feb  4 10:48:00 2021

@author: Daniel
"""


function [Ke, fint, fext, stress,epsilon] = elem4n(ue, ex, ey, ep, mp, eq)
% function [Ke, fint, fext, stress] = elem4n(ue, ex, ey, ep, mp, eq)
% Element routine for 4-node isoparametric elements. 
% Can handle non-linear materials without state variables
% Given element displacement ue, calculate element stiffness Ke, internal
% force vector fint and stress. If body load is given, calculate
% contribution to external load vector fext.
% Function is a modified version of CALFEM's plani4e (incl. plani4s/f)
% 
% Input
% ue    [8x1]   Element displacement [u1,u2,...,u8]'
% ex    [1x4]   Element nodal x-coordinates
% ey    [1x4]   Element nodal y-coordinates
% ep    [1x4]   Element properties: 
%       ep(1):  Analysis type (only 2=plane strain supported)
%               If ep(1)=3 then Abaqus' "Selectively Reduced Integration"
%               will be used (requires ep(3)=2 as well). 
%       ep(2):  Element thickness
%       ep(3):  Integration rule (number of gauss points in each direction)
%       ep(4):  Material model: 1=elasticity, 2=hencky plasticity
% mp    [1xN]   Material parameters for material model
% eq    [2x1]   Body load (per volume) [bx;by] (Optional input)
% 
% Output
% Ke    [8x8]   Element stiffness matrix
% fint  [8x1]   Internal element load vector
% fext  [8x1]   External element load vector
% stress[6xngp] Element stress in each gauss point (1,2,...,ngp)
%               Each column: [sxx; syy; szz; sxy; sxz; syz]
% 
% Written by Knut Andreas Meyer
% FEM Structures 2017
% Modifications
% 2017-01-30:   Changed error in comments
% 2019-03-04:   Added support for selectively reduced integration
% =========================================================================

%% Read input data
ptype   = ep(1);          % Which analysis type?
t       = ep(2);              % Element thickness
ir      = ep(3);  ngp=ir*ir; % Integration rule and number of gauss points
matmod  = ep(4);

% If 6th input argument present, assign body load
if nargin==6,   b=eq;   else   b=zeros(2,1);   end

%% Setup gauss quadrature
[wp, xsi, eta] = gauss_quadrature(ir);

%% Setup shape functions and its derivatives at each gauss point
[N, dNr] = shape_functions(eta, xsi, ngp);
JT=dNr*[ex;ey]';

%% Loop over integration points, add the contributions to Ke, fint and fext
Ke      = zeros(8,8);   %Preallocate Ke
fint    = zeros(8,1);   %Preallocate fint
fext    = zeros(8,1);   %Preallocate fext
stress  = zeros(6, ngp);%Preallocate stress

%Plane strain, selectively reduced integration (Abaqus' method)
bbar = ptype==3 && ir==2;
if bbar
    [~, dNr_bb] = shape_functions(0, 0, 1);
    JT_bb=dNr_bb*[ex;ey]';
    dNx_bb=JT_bb\dNr_bb;
    tmp = zeros(1,8);
    tmp(1:2:end) = dNx_bb(1,:)/3;
    tmp(2:2:end) = dNx_bb(2,:)/3;
    Bvol0 = [tmp; tmp; 0*tmp];
    ptype=2;
end

if ptype==2 %Plane strain
    
    for i=1:ngp
        indx=[ 2*i-1; 2*i ];
        detJ=det(JT(indx,:));
        if detJ<10*eps
            disp('Jacobideterminant equal or less than zero!')
        end
        dNx=JT(indx,:)\dNr(indx,:);
        
%       Extract values of B(xsi, eta) at current gauss point
        B(1,1:2:7) = dNx(1,:);
        B(2,2:2:8) = dNx(2,:);
        B(3,1:2:7) = dNx(2,:);
        B(3,2:2:8) = dNx(1,:);
        
%       If requested, make the bbar correction
        if bbar
            Bvol = zeros(3,8);
            tmp = (B(1,:)+B(2,:))/3;
            Bvol(1:2, :) = [tmp; tmp];
            B = Bvol0 + (B-Bvol);
        end

%       Extract shape function N(xsi, eta) at current gauss point
        N2(1,1:2:7) = N(i,:);
        N2(2,2:2:8) = N(i,:);
        
%       Calculate strain at current gauss point
        epsilon = zeros(6,1);
        epsilon([1;2;4]) = B*ue;
        
%       Calculate material response at current gauss point
        if matmod==1        %Elasticity
            [sigma, dsde] = elastic(epsilon, mp);
        elseif matmod==2    %Hencky plasticity
            [sigma, dsde] = hencky(epsilon, mp);
        elseif matmod==3    %Hencky plasticity
            [sigma, dsde] = Hencky(epsilon, mp);
        else
            error('Only material model (ep(4) 1 or 2 supported');
        end
        
        stress(:, i) = sigma;   %Save stress for current gauss point
        
%       Calculate the gauss point's contribution to element stiffness and forces
        Dm=dsde([1 2 4],[1 2 4]);                  % Components for plane strain
        Ke=Ke+B'*Dm*B*detJ*wp(i)*t;                % Stiffness contribution
        fint=fint+B'*sigma([1;2;4])*wp(i)*detJ*t;  % Internal force vector 
        fext=fext+N2'*b*detJ*wp(i)*t;              % External force vector
    end
else
    error('Only plane strain ep(1)=ptype=2 allowed (unless ep(2)=2, then ep(1)=3 is allowed)');
end

end

function [wp, xsi, eta] = gauss_quadrature(ir)
%% Setup gauss quadrature
if ir==1
    g1=0.0; w1=2.0;
    gp=[ g1 g1 ];  w=[ w1 w1 ];
elseif ir==2
    g1=0.577350269189626; w1=1;
    gp(:,1)=[-g1; g1;-g1; g1];  gp(:,2)=[-g1;-g1; g1; g1];
    w(:,1)=[ w1; w1; w1; w1];   w(:,2)=[ w1; w1; w1; w1];
elseif ir==3
    g1=0.774596669241483; g2=0.;
    w1=0.555555555555555; w2=0.888888888888888;
    gp(:,1)=[-g1;-g2; g1;-g1; g2; g1;-g1; g2; g1];
    gp(:,2)=[-g1;-g1;-g1; g2; g2; g2; g1; g1; g1];
    w(:,1)=[ w1; w2; w1; w1; w2; w1; w1; w2; w1];
    w(:,2)=[ w1; w1; w1; w2; w2; w2; w1; w1; w1];
else
    error('Invalid integration rule. ir = 1, 2 or 3');
end

wp=w(:,1).*w(:,2);
xsi=gp(:,1);  
eta=gp(:,2);

end

function [N, dNr] = shape_functions(eta, xsi, ngp)
    r2=ngp*2;
    N(:,1)=(1-xsi).*(1-eta)/4;  
    N(:,2)=(1+xsi).*(1-eta)/4;
    N(:,3)=(1+xsi).*(1+eta)/4;  
    N(:,4)=(1-xsi).*(1+eta)/4;

    dNr(1:2:r2,1)=-(1-eta)/4;     
    dNr(1:2:r2,2)= (1-eta)/4;
    dNr(1:2:r2,3)= (1+eta)/4;     
    dNr(1:2:r2,4)=-(1+eta)/4;
    
    dNr(2:2:r2+1,1)=-(1-xsi)/4;   
    dNr(2:2:r2+1,2)=-(1+xsi)/4;
    dNr(2:2:r2+1,3)= (1+xsi)/4;   
    dNr(2:2:r2+1,4)= (1-xsi)/4;
    
end