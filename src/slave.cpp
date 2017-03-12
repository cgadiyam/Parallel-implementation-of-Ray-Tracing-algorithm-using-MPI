//This file contains the code that the master process will execute.

#include <iostream>
#include <mpi.h>
#include <math.h>       /* sqrt */
#include "RayTrace.h"
#include "slave.h"

void slaveMain(ConfigData* data)
{
    //Depending on the partitioning scheme, different things will happen.
    //You should have a different function for each of the required 
    //schemes that returns some values that you need to handle.
    switch (data->partitioningMode)
    {
        case PART_MODE_STATIC_STRIPS_HORIZONTAL:
			static_strips_horizontal(data);
            break;
			
		case PART_MODE_STATIC_CYCLES_VERTICAL:
			static_cycles_vertical(data);
            break;
			
        case PART_MODE_STATIC_BLOCKS:
            static_blocks(data);
            break;
			
        case PART_MODE_DYNAMIC:
            dynamic(data);
            break;
			
        default:
            std::cout << "This mode (" << data->partitioningMode;
            std::cout << ") is not currently implemented." << std::endl;
            break;
    }
}

void static_strips_horizontal(ConfigData* data)
{
	float* send_pixels = new float[(3 * data->width * data->height)+1];
    //Start the computation time timer.
    double computationStart = MPI_Wtime();
	
	//calculate rows per strip
	int rowsperblock = (data->height)/(data->mpi_procs);
	int remainingrows = (data->height)%(data->mpi_procs);
	if(data->mpi_rank < remainingrows)
	{
		rowsperblock++;
	}
	
	//assignment of one strip to each slave process
	int offset = 0;
	for(int i = 0;i<data->mpi_rank;i++)
	{
		if(i < remainingrows)
		{
			offset += ((data->height)/(data->mpi_procs)) + 1;
		}
		else
		{
			offset += (data->height)/(data->mpi_procs);
		}
	}
	
    //Render the scene.
    for( int i = 0; i < rowsperblock; ++i)
    {
        for( int j = 0; j < data->width; ++j )
        {
            int row = i+offset;
            int column = j;

            //Calculate the index into the array.
            int baseIndex = 3 * ( row * data->width + column );

            //Call the function to shade the pixel.
            shadePixel(&(send_pixels[baseIndex]),row,j,data);
        }
    }
	
	//Stop the comp. timer
    double computationStop = MPI_Wtime();
    double computationTime = computationStop - computationStart;
	send_pixels[3 * data->width * data->height] = (float)computationTime;
	
	MPI_Send(send_pixels, ((3 * data->width * data->height)+1), MPI_FLOAT, 0, data->mpi_rank, MPI_COMM_WORLD);
}

void static_cycles_vertical(ConfigData* data)
{
	float* send_pixels = new float[(3 * data->width * data->height)+1];
    //Start the computation time timer.
    double computationStart = MPI_Wtime();
	
    //Render the scene.
    for( int j = (data->mpi_rank * data->cycleSize); j < data->width; j += (data->mpi_procs * data->cycleSize))
    {
		for( int k = 0; k < data->cycleSize; ++k )
		{
			if(j+k < data->width)
			{
				for( int i = 0; i < data->height; ++i )
				{
					int row = i;
					int column = j+k;

					//Calculate the index into the array.
					int baseIndex = 3 * ( row * data->width + column );

					//Call the function to shade the pixel.
					shadePixel(&(send_pixels[baseIndex]),row,column,data);
				}
			}
		}
    }
	
	//Stop the comp. timer
    double computationStop = MPI_Wtime();
    double computationTime = computationStop - computationStart;
	send_pixels[3 * data->width * data->height] = (float)computationTime;
	
	MPI_Send(send_pixels, ((3 * data->width * data->height)+1), MPI_FLOAT, 0, data->mpi_rank, MPI_COMM_WORLD);
}

void static_blocks(ConfigData* data)
{
	float* send_pixels = new float[(3 * data->width * data->height)+1];
    //Start the computation time timer.
    double computationStart = MPI_Wtime();
	
	//calculate total number of square blocks and find dimensions
	int temp = (data->width * data->height)/data->mpi_procs;
	int block_length = (int)sqrt(temp);
	int xblocks, yblocks;
	xblocks = yblocks = sqrt(data->mpi_procs);
	int row_offset = 0;
	int column_offset = 0;
	
	//assign one block to each slave process
	for(int i = 0;i<data->mpi_rank;i++)
	{
		column_offset += block_length;
		if(column_offset > (xblocks * block_length)-1)
		{
			column_offset = 0;
			row_offset += block_length;
		}
	}
	
    //Render the scene.
    for( int j = 0; j < block_length; ++j )
    {
        for( int i = 0; i < block_length; ++i )
        {
            int row = i + row_offset;
            int column = j + column_offset;

            //Calculate the index into the array.
            int baseIndex = 3 * ( row * data->width + column );

            //Call the function to shade the pixel.
            shadePixel(&(send_pixels[baseIndex]),row,column,data);
        }
    }

    //Stop the comp. timer
    double computationStop = MPI_Wtime();
    double computationTime = computationStop - computationStart;
	send_pixels[3 * data->width * data->height] = (float)computationTime;
	
	MPI_Send(send_pixels, ((3 * data->width * data->height)+1), MPI_FLOAT, 0, data->mpi_rank, MPI_COMM_WORLD);
}

void dynamic(ConfigData* data)
{
	MPI_Status status;
	float* send_pixels = new float[(3 * data->width * data->height)+1];
	int slave_block_request = -1;
	int row_offset, column_offset;
	row_offset = 0;
	column_offset = 0;
	double computationStart, computationStop, computationTime;
	int slave_offset[3];
	bool termination = false;
	
	//calculate row offset and column offset of initial block
	for(int i = 1;i<data->mpi_rank;i++)
	{
		column_offset += data->dynamicBlockWidth;
		if(column_offset > (data->width - data->dynamicBlockWidth))
		{
			column_offset = 0;
			row_offset += data->dynamicBlockHeight;
		}
		if(row_offset > (data->height - data->dynamicBlockHeight))
		{
			termination = true;
		}
	}
	
	computationTime = 0;
	while(!termination)
	{
		computationStart = MPI_Wtime();
		//Render the scene.
		for( int j = 0; j < data->dynamicBlockWidth; ++j )
		{
			for( int i = 0; i < data->dynamicBlockHeight; ++i )
			{
				int row = i + row_offset;
				int column = j + column_offset;

				//Calculate the index into the array.
				int baseIndex = 3 * ( row * data->width + column );
				//Call the function to shade the pixel.
				shadePixel(&(send_pixels[baseIndex]),row,column,data);
			}
		}
		computationStop = MPI_Wtime();
		computationTime += (computationStop - computationStart);
		slave_block_request = 0;
		//request master for a new block
		MPI_Send(&slave_block_request, 1, MPI_INT, 0, data->mpi_rank, MPI_COMM_WORLD);
		//receive block coordinates from master
		MPI_Recv(slave_offset, 3, MPI_INT, 0, (data->mpi_rank)+500, MPI_COMM_WORLD, &status);
		if(slave_offset[2] == 0)
		{
			column_offset = slave_offset[0];
			row_offset = slave_offset[1];
		}
		else if(slave_offset[2] == -1)
		{
			//if master sends termination notification, send computed data to master
			send_pixels[(3 * data->width * data->height)] = (float)(computationTime);
			MPI_Send(send_pixels, ((3 * data->width * data->height)+1), MPI_FLOAT, 0, (data->mpi_rank)+600, MPI_COMM_WORLD);
			termination = true;
		}
		else
		{
		}
	}
}