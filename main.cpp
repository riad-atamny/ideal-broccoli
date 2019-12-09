#include <iostream>
#include "switch.h"



int main() {
    Switch controller(1, 100,true);

    int inserted = 0;
    while (inserted < 100) {
        Rule rule(1);
        bool res;
        res = controller.ctrlAddRule(rule);
        if (res) {
            inserted++;
        }
    }
}