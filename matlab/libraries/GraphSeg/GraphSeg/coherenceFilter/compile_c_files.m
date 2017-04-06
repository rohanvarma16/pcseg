% This script will compile all the C files of the registration methods

cd('functions')
files=dir('*.c');
for i=1:length(files)
    filename=[files(i).name];
    disp(['compiling : ' filename]);
    mex(filename,'-v');
end
cd('..');


cd('functions2D')
    disp('compiling : CoherenceFilterStep2D');
    mex CoherenceFilterStep2D.c -v;
cd('..');

cd('functions3D')
    disp('compiling : CoherenceFilterStep3D');
    mex CoherenceFilterStep3D.c -v;
    disp('compiling : diffusion_scheme_3D_non_negativity');
    mex diffusion_scheme_3D_non_negativity.c -v;
    disp('compiling : diffusion_scheme_3D_rotation_invariant');
    mex diffusion_scheme_3D_rotation_invariant.c -v;
    disp('compiling : diffusion_scheme_3D_standard');
    mex diffusion_scheme_3D_standard.c -v;
    disp('compiling : EigenVectors3');
    mex EigenVectors3D.c -v;
    disp('compiling : StructureTensor2DiffusionTensor3DWeickert');
    mex StructureTensor2DiffusionTensor3DWeickert.c -v;
cd('..');


