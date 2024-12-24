//---------------------------------------------------------------------------

#ifndef BinFileH
#define BinFileH

struct ImageData
{
    std::vector<uint16_t> data;
    uint32_t width;
    uint32_t height;
    uint32_t samplesPerPixel;

    uint16_t &at(uint32_t row, uint32_t col, uint32_t channel)
    {
        return data[(row * width + col) * samplesPerPixel + channel];
    }

    const uint16_t &at(uint32_t row, uint32_t col, uint32_t channel) const
    {
        return data[(row * width + col) * samplesPerPixel + channel];
    }
};

ImageData binaryToArray(const std::string &filename);
//---------------------------------------------------------------------------
#endif
