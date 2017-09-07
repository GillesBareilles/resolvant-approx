if isunix
    if not((exist('chebyProj.mexglx') == 3) | (exist('chebyProj.mexa64') == 3))
        mex GCC='/usr/bin/g++-4.7' CXXFLAGS="\$CXXFLAGS -std=c++11" -I"./include/" -largeArrayDims src/chebyProj.cpp
    end
elseif ispc
    if not((exist('chebyProj.mexw32') == 3) | (exist('chebyProj.mexw64') == 3))
        mex CXXFLAGS="\$CXXFLAGS -std=c++11" -I"./include/" -largeArrayDims src/chebyProj.cpp
    end
elseif ismac
    if (exist('chebyProj.mexmaci') ~= 3)
        mex CXXFLAGS="\$CXXFLAGS -std=c++11" -I"./include/" -largeArrayDims src/chebyProj.cpp
    end
end
