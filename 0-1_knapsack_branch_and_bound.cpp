#include <stdio.h>
#include <algorithm>
using namespace std;

struct Node {
    int level = -1;
    int w = 0;
    int v = 0;
    int upper_bound = 0;
    Node *parent = NULL;
    Node *l_child = NULL;
    Node *r_child = NULL;
};

struct Item {
    int w = 0;
    int v = 0;
    double cp = 0;
};

bool cmp_func(Item a, Item b) {
    if (a.cp > b.cp) {
        return 1;
    } else {
        return 0;
    }
}

int main() {
    int max_v = 0, max_w = 0, num_item = 0, i = 0;  
    scanf("%d %d", &max_w, &num_item);
    
    int max_level = num_item - 1;
    Item *item = new Item[num_item];
    
    for (i = 0; i < num_item; ++i) {
        scanf("%d %d", &item[i].v, &item[i].w);
        item[i].cp = (double)item[i].v / item[i].w;
    }
    // sort by cp
    sort(item, item + num_item, cmp_func);

    // initial node
    Node *curr_node = new Node();
    int l_flag = 0, level = 0, tmp_v = 0, remain_w = 0;

    while (true) {
        level = curr_node -> level; // level is curr_node -> level
        Node *next_node = new Node();
        
        if ((curr_node -> l_child) == NULL && level != max_level) {
            // build left child
            level++; // now level is next_node -> level
            next_node -> level = level;
            next_node -> w = curr_node -> w + item[level].w;
            next_node -> v = curr_node -> v + item[level].v;
            next_node -> parent = curr_node;
            curr_node -> l_child = next_node;
            l_flag = 1;
        } else if ((curr_node -> r_child) == NULL && level != max_level) {
            // build right child
            level++; // now level is next_node -> level
            next_node -> level = level;
            next_node -> w = curr_node -> w;
            next_node -> v = curr_node -> v;
            next_node -> parent = curr_node;
            curr_node -> r_child = next_node;
            l_flag = 0;
        } else {
            if ((curr_node -> parent) != NULL) {
                delete (curr_node -> r_child);
                delete (curr_node -> l_child);
                curr_node = curr_node -> parent;
                continue;
            } else {
                // head node
                break;
            }
        }

        // now level is next_node -> level
        if ((next_node -> w) <= max_w) {
            if (level != max_level) {
                remain_w = max_w - (curr_node -> w);
                if (remain_w == 0) {
                    continue;
                }

                // compute upper bound, use greedy(frac)
                tmp_v = curr_node -> v;
                if (l_flag == 0) {
                    // right child skip item[level]
                    level++;
                } 
                for (i = level; i <= max_level; ++i) {
                    if (remain_w < item[i].w) {
                        // knapsack's size is not enough
                        tmp_v += (remain_w * item[i].cp);  
                        break;
                    } else {
                        tmp_v += item[i].v;
                        remain_w -= item[i].w;
                    }
                }
                next_node -> upper_bound = tmp_v;

                if ((next_node -> v) > max_v) {
                    max_v = next_node -> v;
                }
                // detemine whether need to check next_node
                if ((next_node -> upper_bound) > max_v) {
                    curr_node = next_node;
                }
            } else {
                // bottom level, check value only
                if ((next_node -> v) > max_v) {
                    max_v = next_node -> v;
                }
            }
        }
    } // while

    printf("%d\n", max_v);
    return 0;
}
