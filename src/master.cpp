//This file contains the code that the master process will execute.

#include <iostream>
#include <mpi.h>
#include <math.h>       /* sqrt */
#include "RayTrace.h"
#include "master.h"

void masterMain(ConfigData* data)
{
    //Depending on the partitioning scheme, different things will happen.
    //You should have a different function for each of the required 
    //schemes that returns some values that you need to handle.
    
    //Allocate space for the image on the master.
    float* pixels = new float[3 * data->width * data->height];
    
    //Execution time will be defined as how long it takes
    //for the given function to execute based on partitioning
    //type.
    double renderTime = 0.0, startTime, stopTime;

	//Add the required partitioning methods here in the case statement.
	//You do not need to handle all cases; the default will catch any
	//statements that are not specified. This switch/case statement is the
	//only place that you should be adding code in this function. Make sure
	//that you update the header files with the new functions that will be
	//called.
	//It is suggested that you use the same parameters to your functions as shown
	//in the sequential example below.
    switch (data->partitioningMode)
    {
        case PART_MODE_STATIC_STRIPS_HORIZONTAL:
            startTime = MPI_Wtime();
            static_strips_horizontal(data, pixels);
            stopTime = MPI_Wtime();
            break;
			
        case PART_MODE_STATIC_CYCLES_VERTICAL:
            startTime = MPI_Wtime();
            static_cycles_vertical(data, pixels);
            stopTime = MPI_Wtime();
            break;
			
        case PART_MODE_STATIC_BLOCKS:
            startTime = MPI_Wtime();
            static_blocks(data, pixels);
            stopTime = MPI_Wtime();
            break;

        case PART_MODE_DYNAMIC:
            startTime = MPI_Wtime();
            dynamic(data, pixels);
            stopTime = MPI_Wtime();
            break;				
			
        default:
            std::cout << "This mode (" << data->partitioningMode;
            std::cout << ") is not currently implemented." << std::endl;
            break;
    }

    renderTime = stopTime - startTime;
    std::cout << "Execution Time: " << renderTime << " seconds" << std::endl << std::endl;

    //After this gets done, save the image.
    std::cout << "Image will be save to: ";
    std::string file = generateFileName(data);
    std::cout << file << std::endl;
    savePixels(file, pixels, data);

    //Delete the pixel data.
    delete[] pixels; 
}

void masterSequential(ConfigData* data, float* pixels)
{
    //Start the computation time timer.
    double computationStart = MPI_Wtime();

    //Render the scene.
    for( int i = 0; i < data->height; ++i )
    {
        for( int j = 0; j < data->width; ++j )
        {
            int row = i;
            int column = j;

            //Calculate the index into the array.
            int baseIndex = 3 * ( row * data->width + column );

            //Call the function to shade the pixel.
            shadePixel(&(pixels[baseIndex]),row,j,data);
        }
    }

    //Stop the comp. timer
    double computationStop = MPI_Wtime();
    double computationTime = computationStop - computationStart;

    //After receiving from all processes, the communication time will
    //be obtained.
    double communicationTime = 0.0;

    //Print the times and the c-to-c ratio
	//This section of printing, IN THIS ORDER, needs to be included in all of the
	//functions that you write at the end of the function.
    std::cout << "Total Computation Time: " << computationTime << " seconds" << std::endl;
    std::cout << "Total Communication Time: " << communicationTime << " seconds" << std::endl;
    double c2cRatio = communicationTime / computationTime;
    std::cout << "C-to-C Ratio: " << c2cRatio << std::endl;
}

void static_strips_horizontal(ConfigData* data, float* pixels)
{
	MPI_Status status;
	float* recv_pixels = new float[(3 * data->width * data->height)+1];
    //Start the computation time timer.
    double computationStart = MPI_Wtime();
	
	//calculate rows per strip
	int rowsperblock = (data->height)/(data->mpi_procs);
	int remainingrows = (data->height)%(data->mpi_procs);
	if(data->mpi_rank < remainingrows)
	{
		rowsperblock++;
	}
	
	//assignment of strip
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
            shadePixel(&(pixels[baseIndex]),row,j,data);
        }
    }

    //Stop the comp. timer
    double computationStop = MPI_Wtime();
    double computationTime = computationStop - computationStart;
	double communicationStart, communicationStop, communicationTime, temp;
	communicationTime = 0;
	
	//receive computed data from all slave processes
	for(int i = 1;i < data->mpi_procs;i++)
	{
		communicationStart = MPI_Wtime();
		MPI_Recv(recv_pixels, ((3 * data->width * data->height)+1), MPI_FLOAT, i, i, MPI_COMM_WORLD, &status);
		communicationStop = MPI_Wtime();
		temp = communicationStop - communicationStart;
		communicationTime+=temp;
		rowsperblock = (data->height)/(data->mpi_procs);
		remainingrows = (data->height)%(data->mpi_procs);
		if(i < remainingrows)
		{
			rowsperblock++;
		}
		offset = 0;
		for(int k = 0;k<i;k++)
		{
			if(k < remainingrows)
			{
				offset += ((data->height)/(data->mpi_procs)) + 1;
			}
			else
			{
				offset += (data->height)/(data->mpi_procs);
			}
		}
		for( int m = 0; m < rowsperblock; ++m)
		{
			for( int n = 0; n < data->width; ++n )
			{
				int row = m + offset;
				int column = n;
				int index = 3 * ( row * data->width + column );
				pixels[index] = recv_pixels[index];
				pixels[index+1] = recv_pixels[index+1];
				pixels[index+2] = recv_pixels[index+2];
			}
		}
		    
		computationTime+=(double)recv_pixels[3 * data->width * data->height];
	}

    //Print the times and the c-to-c ratio
	//This section of printing, IN THIS ORDER, needs to be included in all of the
	//functions that you write at the end of the function.
    std::cout << "Total Computation Time: " << computationTime << " seconds" << std::endl;
    std::cout << "Total Communication Time: " << communicationTime << " seconds" << std::endl;
    double c2cRatio = communicationTime / computationTime;
    std::cout << "C-to-C Ratio: " << c2cRatio << std::endl;
}

void static_cycles_vertical(ConfigData* data, float* pixels)
{
	MPI_Status status;
	float* recv_pixels = new float[(3 * data->width * data->height)+1];
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
					shadePixel(&(pixels[baseIndex]),row,column,data);
				}
			}
		}
    }

    //Stop the comp. timer
    double computationStop = MPI_Wtime();
    double computationTime = computationStop - computationStart;
	double communicationStart, communicationStop, communicationTime, temp;
	communicationTime = 0;
	//receive computed data from all slave processes
	for(int i = 1;i < data->mpi_procs;i++)
	{
		communicationStart = MPI_Wtime();
		MPI_Recv(recv_pixels, ((3 * data->width * data->height)+1), MPI_FLOAT, i, i, MPI_COMM_WORLD, &status);
		communicationStop = MPI_Wtime();
		temp = communicationStop - communicationStart;
		communicationTime+=temp;
		for(  int n = (i * data->cycleSize); n < data->width; n += (data->mpi_procs * data->cycleSize))
		{
			for( int p = 0; p < data->cycleSize; ++p )
			{
				if(n + p < data->width)
				{
					for( int m = 0; m < data->height; ++m )
					{
						int row = m;
						int column = n + p;
						int index = 3 * ( row * data->width + column );
						pixels[index] = recv_pixels[index];
						pixels[index+1] = recv_pixels[index+1];
						pixels[index+2] = recv_pixels[index+2];
					}
				}
			}
		}	    
		computationTime+=(double)recv_pixels[3 * data->width * data->height];
	}

    //Print the times and the c-to-c ratio
	//This section of printing, IN THIS ORDER, needs to be included in all of the
	//functions that you write at the end of the function.
    std::cout << "Total Computation Time: " << computationTime << " seconds" << std::endl;
    std::cout << "Total Communication Time: " << communicationTime << " seconds" << std::endl;
    double c2cRatio = communicationTime / computationTime;
    std::cout << "C-to-C Ratio: " << c2cRatio << std::endl;
}

void static_blocks(ConfigData* data, float* pixels)
{
	MPI_Status status;
	float* recv_pixels = new float[(3 * data->width * data->height)+1];
    //Start the computation time timer.
    double computationStart = MPI_Wtime();
	
	//calculate total number of square blocks and find dimensions
	int temp1 = (data->width * data->height)/data->mpi_procs;
	int block_length = (int)sqrt(temp1);
	int xblocks, yblocks;
	xblocks = yblocks = sqrt(data->mpi_procs);
	int row_offset = 0;
	int column_offset = 0;

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
            shadePixel(&(pixels[baseIndex]),row,j,data);
        }
    }
	
	for( int j = 0; j < data->height; ++j )
    {
        for( int i = (xblocks * block_length); i < data->width; ++i )
        {
            int row = i;
            int column = j;

            //Calculate the index into the array.
            int baseIndex = 3 * ( row * data->width + column );

            //Call the function to shade the pixel.
            shadePixel(&(pixels[baseIndex]),row,j,data);
        }
    }
	
	for( int j = (yblocks * block_length); j < data->height; ++j )
    {
        for( int i = 0; i < (xblocks * block_length); ++i )
        {
            int row = i;
            int column = j;

            //Calculate the index into the array.
            int baseIndex = 3 * ( row * data->width + column );

            //Call the function to shade the pixel.
            shadePixel(&(pixels[baseIndex]),row,column,data);
        }
    }

    //Stop the comp. timer
    double computationStop = MPI_Wtime();
    double computationTime = computationStop - computationStart;
	double communicationStart, communicationStop, communicationTime, temp2;
	communicationTime = 0;
	//receive computed data from all slave processes
	for(int i = 1;i < data->mpi_procs;i++)
	{
		communicationStart = MPI_Wtime();
		MPI_Recv(recv_pixels, ((3 * data->width * data->height)+1), MPI_FLOAT, i, i, MPI_COMM_WORLD, &status);
		communicationStop = MPI_Wtime();
		temp2 = communicationStop - communicationStart;
		communicationTime+=temp2;
		column_offset = row_offset = 0;
		for(int k = 0;k<i;k++)
		{
			column_offset += block_length;
			if(column_offset > (xblocks * block_length)-1)
			{
				column_offset = 0;
				row_offset += block_length;
			}
		}
		for(  int n = 0; n < block_length; ++n )
		{
			for( int m = 0; m < block_length; ++m )
			{
				int row = m + row_offset;
				int column = n + column_offset;
				int index = 3 * ( row * data->width + column );
				pixels[index] = recv_pixels[index];
				pixels[index+1] = recv_pixels[index+1];
				pixels[index+2] = recv_pixels[index+2];
			}
		}	    
		computationTime+=(double)recv_pixels[3 * data->width * data->height];
	}

    //Print the times and the c-to-c ratio
	//This section of printing, IN THIS ORDER, needs to be included in all of the
	//functions that you write at the end of the function.
    std::cout << "Total Computation Time: " << computationTime << " seconds" << std::endl;
    std::cout << "Total Communication Time: " << communicationTime << " seconds" << std::endl;
    double c2cRatio = communicationTime / computationTime;
    std::cout << "C-to-C Ratio: " << c2cRatio << std::endl;
}

void dynamic(ConfigData* data, float* pixels)
{
	float* recv_pixels = new float[(3 * data->width * data->height)+1];
    double computationTime = 0;
	double communicationStart, communicationStop, communicationTime;
	communicationTime = 0;
	int row_offset = 0;
	int column_offset = 0;
	int i = 1;
	int slave_offset[3];
	bool termination = false;
	bool slave_rendering_done = false;
	int term_count = 0;
	int prev_xoffset = 0;
	int prev_yoffset = 0;
	MPI_Status status;
	bool queue_empty = false;
	int active_procs = 0;
	int xval, yval;
	xval = 0;
	yval = 0;
	int num_blocks = 0;
	int slave_block_request = -1;

	//calculate total number of possible blocks in the rendered image and create block list for all slave processes
	num_blocks = ((int)(data->width/data->dynamicBlockWidth)) * ((int)(data->height/data->dynamicBlockHeight));
	int *slave_block_list = new int[(data->mpi_procs - 1) * num_blocks * 2];
	int *array_offset = new int[(data->mpi_procs - 1)];
	for(int i = 0;i<(data->mpi_procs - 1);i++)
	{
		array_offset[i] = 0;
	}
	slave_block_list[0] = prev_xoffset;
	slave_block_list[1] = prev_yoffset;
	array_offset[0]++;
	
	//update slave block list with initial assignment and find next element in the queue
	for(int i = 1;i<(data->mpi_procs - 1);i++)
	{
		if(queue_empty == false)
		{
			prev_xoffset += data->dynamicBlockWidth;
			if(prev_xoffset > (data->width - data->dynamicBlockWidth))
			{
				prev_xoffset = 0;
				prev_yoffset += data->dynamicBlockHeight;
			}
			if(prev_yoffset > (data->height - data->dynamicBlockHeight))
			{
				prev_yoffset = prev_yoffset - data->dynamicBlockHeight;
				if(queue_empty == false)
				{
					slave_rendering_done = true;
					queue_empty = true;
					active_procs = i;
				}
			}
			else
			{
				slave_block_list[(i * num_blocks * 2)] = prev_xoffset;
				slave_block_list[(i * num_blocks * 2) + 1] = prev_yoffset;
				array_offset[i]++;
			}
			
		}
	}

	//assign blocks dynamically until termination and receive computed data from slave processes upon termination
	while(!termination)
    {
		communicationStart = MPI_Wtime();
		MPI_Recv(&slave_block_request, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		communicationStop = MPI_Wtime();
		communicationTime += (communicationStop - communicationStart);
		i = status.MPI_SOURCE;
		if(slave_rendering_done == false)
		{
			//calculate next element in the queue
			prev_xoffset += data->dynamicBlockWidth;
			if(prev_xoffset > (data->width - data->dynamicBlockWidth))
			{
				prev_yoffset += data->dynamicBlockHeight;
				prev_xoffset = 0;
			}
			if(prev_yoffset > (data->height - data->dynamicBlockHeight))
			{
				slave_rendering_done = true;
			}
		}
		if(slave_rendering_done == false)
		{
			//send coordinates of new block
			slave_offset[0] = prev_xoffset;
			slave_offset[1] = prev_yoffset;
			slave_offset[2] = 0;
			communicationStart = MPI_Wtime();
			MPI_Send(slave_offset, 3, MPI_INT, i, i+500, MPI_COMM_WORLD);
			communicationStop = MPI_Wtime();
			communicationTime += (communicationStop - communicationStart);
			slave_block_list[((i-1) * num_blocks * 2) + (array_offset[i-1] * 2)] = prev_xoffset;
			slave_block_list[((i-1) * num_blocks * 2) + (array_offset[i-1] * 2) + 1] = prev_yoffset;
			array_offset[i-1]++;
		}
		else
		{
			//send termination notification and receive computed data
			slave_offset[2] = -1;
			communicationStart = MPI_Wtime();
			MPI_Send(slave_offset, 3, MPI_INT, i, i+500, MPI_COMM_WORLD);
			MPI_Recv(recv_pixels, ((3 * data->width * data->height)+1), MPI_FLOAT, i, i+600, MPI_COMM_WORLD, &status);
			communicationStop = MPI_Wtime();
			communicationTime += (communicationStop - communicationStart);
			term_count++;
			for(int k = 0;k < array_offset[i-1];k++)
			{
				column_offset = slave_block_list[((i-1) * num_blocks * 2) + (k * 2)];
				row_offset = slave_block_list[((i-1) * num_blocks * 2) + (k * 2) + 1];
				for(  int n = 0; n < data->dynamicBlockWidth; ++n )
				{
					for( int m = 0; m < data->dynamicBlockHeight; ++m )
					{
						int row = m + row_offset;
						int column = n + column_offset;
						int index = 3 * ( row * data->width + column );
						pixels[index] = recv_pixels[index];
						pixels[index+1] = recv_pixels[index+1];
						pixels[index+2] = recv_pixels[index+2];
					}
				}
			}
			computationTime+=(double)(recv_pixels[3 * data->width * data->height]);
		}

		if(queue_empty == false)
		{
			if(term_count == (data->mpi_procs - 1))
			{
				termination = true;
			}
		}
		else
		{
			if(term_count == active_procs)
			{
				termination = true;
			}
		}
    }
	
	bool exit = false;
	prev_xoffset = 0;
	while(!exit)
	{
		prev_xoffset += data->dynamicBlockWidth;
		if(prev_xoffset > (data->width - data->dynamicBlockWidth))
		{
			exit = true;
		}
	}
	
	//do leftover work
	double computationStart = MPI_Wtime();
	//column strip at extreme right
	for( int j = 0; j < data->height; ++j )
    {
        for( int i = prev_xoffset; i < data->width; ++i )
        {
            int row = j;
            int column = i;

            //Calculate the index into the array.
            int baseIndex = 3 * ( row * data->width + column );

            //Call the function to shade the pixel.
            shadePixel(&(pixels[baseIndex]),row,column,data);
        }
    }
	//row strip at the top
	for( int j = prev_yoffset; j < data->height; ++j )
    {
        for( int i = 0; i < prev_xoffset; ++i )
        {
            int row = j;
            int column = i;

            //Calculate the index into the array.
            int baseIndex = 3 * ( row * data->width + column );

            //Call the function to shade the pixel.
            shadePixel(&(pixels[baseIndex]),row,column,data);
        }
    }

    //Stop the comp. timer
    double computationStop = MPI_Wtime();
	//std::cout << "\nmaster: "<<" computation time: "<<(double)(computationStop - computationStart)<<" checkpoint 8 ... ";
    computationTime += (computationStop - computationStart);

    //Print the times and the c-to-c ratio
	//This section of printing, IN THIS ORDER, needs to be included in all of the
	//functions that you write at the end of the function.
    std::cout << "Total Computation Time: " << computationTime << " seconds" << std::endl;
    std::cout << "Total Communication Time: " << communicationTime << " seconds" << std::endl;
    double c2cRatio = communicationTime / computationTime;
    std::cout << "C-to-C Ratio: " << c2cRatio << std::endl;
}