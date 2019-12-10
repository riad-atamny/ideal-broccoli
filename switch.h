//
// Created by ASUS on 24/11/2019.
//

#ifndef COOPERATIVERULECACHING_SWITCH_H
#define COOPERATIVERULECACHING_SWITCH_H

#include <list>
#include "rule.h"

#define MAXRULES 10000

//#include "Controller.h"
typedef enum {

}Actions;
using namespace std;
class Switch {
protected:
    int id;
    int cache_size;
    list<Rule> rules;

public:
    Switch(int id, int cache_size):id(id),cache_size(cache_size){};
    virtual ~Switch()= default;

    virtual Rule help(Packet pkt);
    virtual bool AddRule(Rule rule);

//    Actions action(Packet packet);

    bool proccessData(Packet packet, Switch* controller) {
        // do i have a rule for the packet
        // if yes Salamtak
        for (auto it=rules.begin();it != rules.end();it++){
            if((*it).packetMatch(packet)){
                //do action; for now it is just true!
                return true;
            }
        }
        // if no go to comrade controller
        Rule newRule = controller->help(packet);
        controller->AddRule(newRule);
        // if there is room - Put into cache
        // if not
        //remove least popular rule
        // insert instead

    }

};
Rule Switch::help(Packet pkt) {
    /*run over the whole list of rules till you find the matching rule
     * **this is for cooperative** figure the switch_id and take it from there
     *non-cooperative
     * eliminate the least cached rule and put the new one instead
     * */
    for(auto it=rules.begin();it!=rules.end();it++){
        if((*it).packetMatch(pkt)){
            return *it;
        }
    }
}


#endif //COOPERATIVERULECACHING_SWITCH_H
