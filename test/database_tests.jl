using DFControl

@testset "environment" begin
    e = Environment("environment1", "mpirun -np 2", ["#SBATCH --partition=default", "#SBATCH --time=23:59:00"], ["OMP_NUM_THREADS=1"])
    save(test_server, e)
    e1 = load(test_server, Environment("environment1"))
    @test e == e1
    e1.name = ""

    e1 = load(test_server, e1)
    @test e1 == e
    e1.MPI_command = "mpirun -np 5"
    save(test_server, e1)
    @test load(test_server, Environment("environment1")).MPI_command == "mpirun -np 5"
    e1.name = "environment2"
    e1.MPI_command = "mpirun -np 6"
    save(test_server, e1)
    e1.name = "environment3"
    e1.exports = String[]
    save(test_server, e1)
    @test length(load(test_server, Environment(MPI_command = "mpirun -np 6"))) == 2
    @test length(load(test_server, Environment(MPI_command = "mpirun -np 6", exports=["OMP_NUM_THREADS=1"]))) == 1

    rm(test_server, Environment("environment1"))
    rm(test_server, Environment("environment2"))
    rm(test_server, Environment("environment3"))
    @test length(load(test_server, Environment(""))) == 1
end
