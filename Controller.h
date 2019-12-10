//
// Created by ASUS on 09/12/2019.
//

#ifndef COOPERATIVERULECACHING_CONTROLLER_H
#define COOPERATIVERULECACHING_CONTROLLER_H

#include "switch.h"

class Controller :public Switch{
    list<Switch> switches;

public:
    Controller(int id1, int cacheSize, bool ctrl, int id, int cache_size) : Switch(id1, cacheSize) {
        id=id;

    };
    Packet sendTraffic(int switch_id) {
        auto it = rules.begin();
        int sentRuleI = rand() % MAXRULES;
        for (int i=0; i < sentRuleI; i++) {
            it++;
        }
        bitset<128> pkt =(*it).pattern;
        return Packet(pkt);
    }

    Rule help(Packet pkt) override ;
    bool AddRule(Rule rule){
        if(!(rules.size()<cache_size)) {
            float min=1;
            int i=0;
            _List_iterator<Rule> localRule=rules.begin(),iterator;
            for (;localRule!=rules.end();localRule++) {
                //TODO: should be changed for an array of lambdas in cooperative blah blah!
                if((*localRule).lambda[0]<min){
                    min=(*localRule).lambda[0];
                    iterator=localRule;
                }
            }
            rules.remove(*iterator);
        }

        for (Rule localRule : rules) {
            if (localRule.pattern == rule.pattern) {
                return false;
            }
        }
        rules.push_back(rule);
        rules.sort();
        return true;
    }
};
#endif //COOPERATIVERULECACHING_CONTROLLER_H
