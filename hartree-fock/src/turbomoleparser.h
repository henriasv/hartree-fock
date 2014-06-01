#ifndef TURBOMOLEPARSER_H
#define TURBOMOLEPARSER_H

#include <primitive.h>
#include <contracted.h>
#include <vector>
#include <fstream>
#include <iostream>
#include <armadillo>

class TurboMoleParser
{
public:
    TurboMoleParser();
    TurboMoleParser(char* filename);
    ~TurboMoleParser();
    std::vector<Contracted> returnContracted(arma::vec nucleusPosition);
    void setFile(char* filename);
private:
    void createContracted(arma::vec nucleusPosition);
    std::vector<Contracted> m_contracted;
    std::ifstream m_inFile;
    bool startswith(char[], char[]);
};

#endif // TURBOMOLEPARSER_H
