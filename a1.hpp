/*  KRISHNA PRASAD
 *  PORANDLA
 *  KPORANDL
 */

#ifndef A1_HPP
#define A1_HPP

#include <vector>


void isort(std::vector<short int>& Xi, MPI_Comm comm) {

    int size, rank;

    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);

    //local_counts to store the counts elements assigned in each processor
    std::vector<long long int> local_counts(2*size,0);
    //reduction_sum_counts to store the counts of elements during reduction
    std::vector<long long int> reduction_sum_counts(2*size,0);
    //final_counts_per_processor to store counts of elements that need to present on respective processor
    std::vector<long long int> final_counts_per_processor(2, 0);

    //Performing local count of elements O(n/p)
    //(if total p=4, In each processor local_counts[0] will have count of -3, local_counts[1] will have count of -2 and so on...)
    for(long long int i=0;i<Xi.size();i++){
        local_counts[(size-1)+Xi[i]] += 1;
    }

    //Reduce the local_counts and store in last processor
    MPI_Reduce(local_counts.data(),reduction_sum_counts.data(), local_counts.size(), MPI_LONG_LONG_INT, MPI_SUM, size-1, comm);

    //Scatter back the summed-up counts to each processor
    MPI_Scatter(reduction_sum_counts.data(),2,MPI_LONG_LONG_INT, final_counts_per_processor.data(),2,MPI_LONG_LONG_INT, size-1, comm);

    long long int firstElementCount = final_counts_per_processor[0];
    long long int secondElementCount = final_counts_per_processor[1];

    //Xi will be resized as per the total number of elements assigned to current processor
    Xi.resize(firstElementCount+secondElementCount);

    //Calculating first and second elements based on rank and size
    short int firstElement = (-(size-1)) + (2*rank);
    short int secondElement  = firstElement+1;

    //Finally, assigning all the elements to Xi
    for(long long int i=0;i<firstElementCount;i++){
        Xi[i] = firstElement;
    }
    for(long long int i=0;i<secondElementCount;i++){
        Xi[firstElementCount + i] = secondElement;
    }

} // isort

#endif // A1_HPP


    