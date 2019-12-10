#include <iostream>
#include "Controller.h"



int main() {
    Controller controller(0, 0, false, 1, 100);

    int inserted = 0;
    while (inserted < 100) {
        Rule rule(1);
        bool res;
        res = controller.AddRule(rule);
        if (res) {
            inserted++;
        }
    }
}