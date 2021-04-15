@testset "Wannier band Python to Julia Equivalence " begin
    kpoints = rand(1, 3)
    reshape_test_array = rand(100)
    reshape_test_array1 = rand(2, 10, 10)
    reshape_test_array2 = rand(4, 2, 8, 9, 10)
    conjugate_test_array = rand(ComplexF64, 5, 5)
    np = pyimport("numpy")
    cell_map =rand(10, 3)
    @test cell_map*(kpoints[1, :]) ≈ np.ndarray.flatten(np.dot(kpoints, np.transpose(cell_map)))
    #@test np.tensordot(np.exp((2j*np.pi)*np.dot(kpoints, np.transpose(cellMap))), Hwannier, axes=1)
    @test np.pi ≈ pi

    @test_broken reshape(reshape_test_array, (25, 4)) ≈ np.reshape(reshape_test_array, (25, 4))
    @test permutedims(reshape(reshape_test_array, (10, 10)), (2, 1)) ≈ np.reshape(reshape_test_array, (10, 10))
    @test permutedims(reshape(np.reshape(reshape_test_array1, (2, 10*10)), (2, 10, 10)), [1, 3, 2]) ≈ reshape_test_array1
    @test np.exp(np.array([1, 2, 3, 4])) ≈ exp.([1, 2, 3, 4])
    @test np.conj(conjugate_test_array) ≈ conj(conjugate_test_array) #Test that the conjugation method in Julia is equivalent to that in python
    @test vec(permutedims(reshape_test_array2, (5, 4, 3, 2, 1))) == np.ndarray.flatten(reshape_test_array2)
end
