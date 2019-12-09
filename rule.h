//
// Created by ASUS on 24/11/2019.
//

#ifndef COOPERATIVERULECACHING_RULE_H
#define COOPERATIVERULECACHING_RULE_H

#include "packet.h"
#include <bitset>

std::bitset<128> generatePattern() {
    std::bitset<128> ans;
    for (int i = 0 ; i < ans.size(); i++) {
        ans[i] = rand() % 2;
    }
    return ans;
}

using namespace std;
class Rule {
public:
    int switch_id;
    float lambda[1];
    int dc_num;
    int action;
    int times;
    bitset<128> pattern;

    Rule(int switch_id):switch_id(switch_id){
        pattern=generatePattern();
        dc_num=rand()%129;
        //TODO:action generate
    }
    Rule(int switch_id, bitset<128> pattern, int dc_num);
    ~Rule()= default;
    bool packetMatch(Packet pkt){

    }
    bool operator<(Rule rule){
        // it is meant to be sorted in descending sort by dc_num
        if( this->dc_num>rule.dc_num){
            return true;
        } else if(this->dc_num==rule.dc_num){
            return lambda[this->switch_id]>lambda[rule.switch_id];
        }
        return false;
    }




};
#endif //COOPERATIVERULECACHING_RULE_H
