function [Chrom,ObjV] = rog_kx_ky(ROGR,NIND,NVAR,PRECI,ObjV,FieldD,Chrom,ObjfunArg)
% ROG   random offspring generation
%         offsprings are inserted in place of weakest individuals
%         for objective function in kx_ky plane
% 
% Syntax: [Chrom,ObjV] = rog(ROGR,NIND,NVAR,PRECI,ObjV,FieldD,Chrom,ObjfunArg) 
% 
% Inputs: 
%    ROGR - random offspring generation rate, double 
%    NIND - number of individuals, integer
%    NVAR - number of variables in objective function, integer
%    PRECI - Precision of binary representation of variables
%    ObjV - objective function value; vector [NIND,1]
%    FieldD - field descriptor - refer to Sheffield GA toolbox
%    Chrom - chromosome to be updated - vector of binary representation of the individuals
%    ObjfunArg - structure containing arguments of objective function
% 
% Outputs: 
%    Chrom - Chromosome with randomly generated offsprings
%    output2 - Description, double, dimensions [m, n], Units: m/s^2 
% 
% Example: 
%    [Chrom,ObjV] = rog(ROGR,NIND,NVAR,PRECI,ObjV,FieldD,Chrom,ObjfunArg) 
%
% Other m-files required: crtbp, bs2rv, obj_ga_C_tensor_known_mass_unidirectional 
% Subfunctions: none 
% MAT-files required: none 
% See also: Sheffield GA Toolbox
% 

% Author: Pawel Kudela, D.Sc., Ph.D., Eng. 
% Institute of Fluid Flow Machinery Polish Academy of Sciences 
% Mechanics of Intelligent Structures Department 
% email address: pk@imp.gda.pl 
% Website: https://www.imp.gda.pl/en/research-centres/o4/o4z1/people/ 

%---------------------- BEGIN CODE---------------------- 

NRANDIND = ceil(ROGR*NIND);
RandChrom = crtbp(NRANDIND, NVAR*PRECI); % NRANDIND random individuals
% evaluate random individual
[ObjVRand] = obj_ga_C_tensor_known_mass_unidirectional_kx_ky(bs2rv(RandChrom,FieldD),ObjfunArg.Data_polar,ObjfunArg.layup,ObjfunArg.h,ObjfunArg.fmin,ObjfunArg.fmax,ObjfunArg.wavenumber_max,ObjfunArg.number_of_frequency_points,ObjfunArg.beta,ObjfunArg.stack_dir,ObjfunArg.np,ObjfunArg.nele_layer,ObjfunArg.number_of_modes_considered,ObjfunArg.rho);
       
[A,I]=sort(ObjV,'descend');
% replace weakest individuals
for k=1:NRANDIND
    Chrom(I(k),:) =  RandChrom(k,:);
    % update ObjV
    ObjV(I(k))=ObjVRand(k);
end


%---------------------- END OF CODE---------------------- 

% ================ [rog.m] ================  
