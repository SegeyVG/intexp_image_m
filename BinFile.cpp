//---------------------------------------------------------------------------

#pragma hdrstop

#include <fstream>
#include <vector>
#include <cstdint>
#include <iostream>
#include <stdexcept>

#include "BinFile.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)

struct BinaryHeader
{
    uint32_t width;           // ������ �����������
    uint32_t height;          // ������ �����������
    uint16_t samplesPerPixel; // ���������� �������
    uint16_t bitsPerSample;   // ��� �� �����
    uint64_t dataOffset;      // �������� ������ �� ������ �����
};




ImageData binaryToArray(const std::string &filename)
{
    std::ifstream inFile(filename, std::ios::binary);
    if (!inFile)
    {
        throw std::runtime_error("������: �� ������� ������� �������� ����: " + filename);
    }

    BinaryHeader header;
    inFile.read(reinterpret_cast<char *>(&header), sizeof(BinaryHeader));
    if (!inFile)
    {
        throw std::runtime_error("������: �� ������� ��������� ���������!");
    }

    if (header.bitsPerSample != 16)
    {
        throw std::runtime_error("������: ����������� �� 16-������!");
    }

    ImageData imgData;
    imgData.width = header.width;
    imgData.height = header.height;
    imgData.samplesPerPixel = header.samplesPerPixel;

    // alloc
    imgData.data.resize(header.width * header.height * header.samplesPerPixel);

    inFile.read(reinterpret_cast<char *>(imgData.data.data()),
                imgData.data.size() * sizeof(uint16_t));
    if (!inFile)
    {
        throw std::runtime_error("������: �� ������� ��������� ������ �����������!");
    }

    inFile.close();
    return imgData;
}