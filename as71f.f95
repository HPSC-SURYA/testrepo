	program AS71f
		include 'mpif.h'
		integer :: my_pe, n_pes, status
		call MPI_Init( status )
		call MPI_Comm_rank( MPI_COMM_WORLD, my_pe, status )
		call MPI_Comm_size( MPI_COMM_WORLD, n_pes, status )
		print *, 'Hello, world from PE ',my_pe,' of ',n_pes
		call MPI_Barrier( MPI_COMM_WORLD, status )
		call MPI_Finalize( status )
	end program AS71f
