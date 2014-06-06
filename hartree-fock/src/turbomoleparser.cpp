#include "turbomoleparser.h"


TurboMoleParser::TurboMoleParser()
{
}

TurboMoleParser::TurboMoleParser(char *filename)
{
    setFile(filename);
}



TurboMoleParser::~TurboMoleParser()
{
    m_inFile.close();
}

std::vector<Contracted> TurboMoleParser::returnContracted(arma::vec nucleusPosition)
{
    createContracted(nucleusPosition);
    return m_contracted;
}

void TurboMoleParser::setFile(char* filename)
{
    m_inFile.open(filename);
    if (!m_inFile) {
        std::cout << "Could not read file "<< filename << std::endl;
        exit(1);
    }
}

void TurboMoleParser::createContracted(arma::vec nucleusPosition)
{
    m_inFile.clear() ;
    m_inFile.seekg(0, std::ios::beg) ;
    m_contracted.clear();
    char buf[256]; // Random size buffer, assuming no important lines to be longer than this
    std::stringstream s;
    int numBasisFunctionsInContracted;
    char orbitalType;
    double coefficient;
    double exponent;

    while (m_inFile.getline(buf, 256)) {
        if (startswith(buf, "$basis")){
            m_inFile.getline(buf, 256);
            m_inFile.getline(buf, 256);
            m_inFile.getline(buf, 256);

            // From here is the information on the basis
            m_inFile.getline(buf, 256);
            while (!startswith(buf, "*")) {
                s << buf;
                s >> numBasisFunctionsInContracted;
                s >> orbitalType;

                switch (orbitalType) {
                case 's':
                {
                    Contracted contracted;
                    for (int numPrimitives = 0; numPrimitives<numBasisFunctionsInContracted; numPrimitives++) {
                        m_inFile.getline(buf, 256);
                        s << buf;
                        s >> exponent;
                        s >> coefficient;
                        Primitive primitive(coefficient, 0, 0, 0, exponent, nucleusPosition);
                        contracted.addPrimitive(primitive);
                    }
                    m_contracted.push_back(contracted);
                    m_inFile.getline(buf, 256);
                    break;
                }
                case 'p':
                {
                    Contracted contractedX; Contracted contractedY; Contracted contractedZ;
                    for (int numPrimitives = 0; numPrimitives<numBasisFunctionsInContracted; numPrimitives++) {
                        m_inFile.getline(buf, 256);
                        s << buf;
                        s >> exponent;
                        s >> coefficient;
                        Primitive px(coefficient, 1, 0, 0, exponent, nucleusPosition);
                        Primitive py(coefficient, 0, 1, 0, exponent, nucleusPosition);
                        Primitive pz(coefficient, 0, 0, 1, exponent, nucleusPosition);
                        contractedX.addPrimitive(px);
                        contractedY.addPrimitive(py);
                        contractedZ.addPrimitive(pz);
                    }
                    m_contracted.push_back(contractedX);
                    m_contracted.push_back(contractedY);
                    m_contracted.push_back(contractedZ);
                    m_inFile.getline(buf, 256);
                    break;
                }
                }
            }
        }
    }
}

bool TurboMoleParser::startswith(char str1[], char str2[])
{
    std::string tmp1 (str1);
    std::string tmp2 (str2);
    return (tmp1.find(str2) == 0);
}

