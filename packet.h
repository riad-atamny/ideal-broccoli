//
// Created by ASUS on 24/11/2019.
//

#ifndef COOPERATIVERULECACHING_PACKET_H
#define COOPERATIVERULECACHING_PACKET_H

#include <bitset>

class Packet{
    std::bitset<128> bits;
public:
    Packet(std::bitset<128> bits){
        this->bits=bits;
    }
    std::bitset<128> packetGetBits(){
        return bits;
    }
};
#endif //COOPERATIVERULECACHING_PACKET_H
