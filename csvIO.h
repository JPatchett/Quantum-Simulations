#ifndef CSVIO_H
#define CSVIO_H

#include <string>
#include <vector>

class CSVWriter{
    
    private:
        std::string filename;
        std::string delimeter;
        int linesCount;
    
    public: 
        CSVWriter(std::string filename, std::string delm = ",");
        void addDataInRow(std::vector<double> & first);
};

#endif