if verLessThan('matlab','9.3')
    mex_command = 'mex';
else
    mex_command = 'mex -R2018a -v';
end

if ispc
    % Windows
   cc = mex.getCompilerConfigurations('C','Selected');
   if(strcmp(cc.Manufacturer,'GNU'))
     mex_first_command = strjoin([mex_command, "CC=gcc"]);
     mex_last_command = 'CFLAGS="$CFLAGS -O3 -mtune=native -mfpmath=both -fopenmp" LDFLAGS="$LDFLAGS -fopenmp -static-libgcc -lgomp"';
     eval(strjoin([mex_first_command, "c_charge.c", mex_last_command]));
     eval(strjoin([mex_first_command, "c_current.c", mex_last_command]));
     eval(strjoin([mex_first_command, "c_rvelocity.c", mex_last_command]));
     eval(strjoin([mex_first_command, "c_energy.c", mex_last_command]));
     eval(strjoin([mex_first_command, "c_calcEJ.c", mex_last_command]));
   else
    eval(strjoin([mex_command, "c_charge.c"]));
    eval(strjoin([mex_command, "c_current.c"]));
    eval(strjoin([mex_command, "c_rvelocity.c"]));
    eval(strjoin([mex_command, "c_energy.c"]));
    eval(strjoin([mex_command, "c_calcEJ.c"]));
   end
    
elseif ismac
    % Windows Mac: Xcode compilers does not support openmp
    eval(strjoin([mex_command, "c_charge.c"]));
    eval(strjoin([mex_command, "c_current.c"]));
    eval(strjoin([mex_command, "c_rvelocity.c"]));
    eval(strjoin([mex_command, "c_energy.c"]));
    eval(strjoin([mex_command, "c_calcEJ.c"]));
elseif isunix
    %Linux: gcc accepts openmp
    mex_first_command = strjoin([mex_command, "CC=gcc"]);
    mex_last_command = 'CFLAGS="$CFLAGS -O3 -mtune=native -fopenmp" LDFLAGS="$LDFLAGS -fopenmp -liomp5"';
    eval(strjoin([mex_first_command, "c_charge.c", mex_last_command]));
    eval(strjoin([mex_first_command, "c_current.c", mex_last_command]));
    eval(strjoin([mex_first_command, "c_rvelocity.c", mex_last_command]));
    eval(strjoin([mex_first_command, "c_energy.c", mex_last_command]));
    eval(strjoin([mex_first_command, "c_calcEJ.c", mex_last_command]));
else
    disp('Platform not supported')
end
