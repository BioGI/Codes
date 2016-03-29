for i in `cat radiusStraightTubeRuns`
do
    sed "s/radiusStraightTube/${i}/" inputSample.txt > input.txt
    mpirun -np 4 Intestine3D.exe &> output
    mkdir multiGridTests/scalarAbsorbedTest/dualLattice/gridRatio4/${i}
    cp scalar*dat multiGridTests/scalarAbsorbedTest/dualLattice/gridRatio4/${i}
done

