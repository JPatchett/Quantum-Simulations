#include "csvIO.h"
#include <fstream>
#include <vector>


CSVWriter::CSVWriter(std::string filename, std::string delm): filename(filename), delimeter(delm), linesCount(0){}

 void CSVWriter::addDataInRow(std::vector<double> & first){

    std::fstream file;

    //Open the file, if linesCount is non-zero then we will add to the existing data, else we will delete it and start a new file
    file.open(filename, std::ios::out| (linesCount == 0? std::ios::trunc : std::ios::app));

    for(int i=0; i< first.size()-1; i++){ //Write data to file;

        file << first[i];
        file << delimeter;
    }

    file << first[first.size()-1];
    file << "\n";
    linesCount++;

    file.close();
}
