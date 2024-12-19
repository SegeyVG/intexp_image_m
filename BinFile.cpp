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
    uint32_t width;           // ширина изображения
    uint32_t height;          // высота изображения
    uint16_t samplesPerPixel; // количество каналов
    uint16_t bitsPerSample;   // бит на канал
    uint64_t dataOffset;      // смещение данных от начала файла
};




ImageData binaryToArray(const std::string &filename)
{
    std::ifstream inFile(filename, std::ios::binary);
    if (!inFile)
    {
        throw std::runtime_error("Ошибка: Не удалось открыть бинарный файл: " + filename);
    }

    BinaryHeader header;
    inFile.read(reinterpret_cast<char *>(&header), sizeof(BinaryHeader));
    if (!inFile)
    {
        throw std::runtime_error("Ошибка: Не удалось прочитать заголовок!");
    }

    if (header.bitsPerSample != 16)
    {
        throw std::runtime_error("Ошибка: Изображение не 16-битное!");
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
        throw std::runtime_error("Ошибка: Не удалось прочитать данные изображения!");
    }

    inFile.close();
    return imgData;
}