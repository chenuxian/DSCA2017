#include <stdio.h>
#include <math.h>
#include <limits.h>

struct linear {
    int x;
    int y;
    int c;
    int next;
};

int rounding(float temp) {
	if (temp > 0) {
		return (int)(temp + 0.5);
	} else if (temp < 0) { 
		return (int)(temp - 0.5);
	} else {
        return 0;
    }
}

int intersect(struct linear* p1, struct linear* p2, float* tmp) {
	// crammer's rule
    int delta = p1->x * p2->y - p2->x * p1->y;
    int delta_x = p1->c * p2->y - p2->c * p1->y;
    int delta_y = p1->x * p2->c - p2->x * p1->c;

    if (delta != 0) {
        (*tmp) = (float)delta_x / delta;
        return 1;
    } else if (delta_x == 0 && delta_y == 0) {
        return 2;
	} else {
        return 3;
	}
}

int intersect2(struct linear* p1, struct linear* p2, float* tmp, float* tmp2) {
    // crammer rule
	int delta = p1->x * p2->y - p2->x * p1->y;
    int delta_x = p1->c * p2->y - p2->c * p1->y;
    int delta_y = p1->x * p2->c - p2->x * p1->c;

    if (delta != 0) {
        (*tmp) = (float)delta_x / delta;
        (*tmp2) = (float)delta_y / delta;
        return 1;
    } else {
        return 0;
	}
}

void override(struct linear* p1, struct linear* p2) {
    p1 -> x = p2 -> x;
    p1 -> y = p2 -> y;
    p1 -> c = p2 -> c;
    p1 -> next = p2 -> next;
}

int main() {
    int i, j, k, line, pos = -1, neg = -1, result;
    scanf("%d", &line);
    struct linear data[line];
    float ax[line];
    float tmp;

    for(i = 0; i < line; ++i) {
        scanf("%d", &data[i].x);
        scanf("%d", &data[i].y);
        scanf("%d", &data[i].c);
        data[i].next = -1;
    }

    // split the data into positive and negative
    for (i = 0; i < line; ++i) {
        if (data[i].y < 0) {
            if (neg < 0) {
                neg = i;
			} else {
                data[k].next = i;
			}
            k = i;
        } else {
            if(pos < 0) {
                pos = i;
			} else {
                data[j].next = i;
			}
            j = i;
        }
    }
    // start Megiddo Algorithm, prune the line into 3
    float xl = -INFINITY, xr = INFINITY;
    while(3 < line) {
        // remove the unnecessary line
        j = -1;
        i = neg;
        while(-1 < i && -1 < data[i].next) {
            --line;
			result = intersect(&data[i], &data[data[i].next], &tmp);
			if (result == 1) {
				// unique intersection
                if(xr < tmp)
                {
					// remove the bigger slope 
                    if((double)data[i].x / data[i].y > (double)data[data[i].next].x / data[data[i].next].y) {
                        data[i].next = data[data[i].next].next;
                    } else {
                        override(&data[i], &data[data[i].next]);
                    }
                }
                else if(tmp < xl)
                {
					// remove the smaller slope
                    if((double)data[i].x / data[i].y > (double)data[data[i].next].x / data[data[i].next].y) {
                        override(&data[i], &data[data[i].next]);
                    } else {
                        data[i].next = data[data[i].next].next;
                    }
                }
                else
                {
					
                    ++line;
                    ax[++j] = tmp;
                    i = data[i].next;
                }
            } else if (result == 2) {
				// infinite intersection
                data[i].next = data[data[i].next].next;
            } else {
				// no intersection
                if(data[i].c > data[data[i].next].c) {
                    data[i].next = data[data[i].next].next;
				} else {
                    override(&data[i], &data[data[i].next]);
				}
            }
            i = data[i].next;
        }
        i = pos;
        while(-1 < i && -1 < data[i].next) {
            --line;
			result = intersect(&data[i], &data[data[i].next], &tmp);
            if (result == 1) {
                if(xr < tmp) {
                    if((double)data[i].x / data[i].y < (double)data[data[i].next].x / data[data[i].next].y) {
                        data[i].next = data[data[i].next].next;
                    } else {
                        override(&data[i], &data[data[i].next]);
                    }
                } else if(tmp < xl) {
                    if((double)data[i].x / data[i].y < (double)data[data[i].next].x / data[data[i].next].y) {
                        override(&data[i], &data[data[i].next]);
                    } else {
                        data[i].next = data[data[i].next].next;
                    }
                } else {
                    ++line;
                    ax[++j] = tmp;
                    i = data[i].next;
                }
            } else if (result == 2) {
                data[i].next = data[data[i].next].next;
			} else {
                if (data[i].c < data[data[i].next].c) {
                    data[i].next = data[data[i].next].next;
				} else {
                    override(&data[i], &data[data[i].next]);
				}
            }
            i = data[i].next;
        }
		
        // find the median of ax, use bubble sort
		int p,q;
		float swap;
		for (p = 0; p < j - 1; ++p) {
    		for (q = 0; q < j - p - 1; q++) {
      			if (ax[q] > ax[q+1]) {
        			swap       = ax[q];
        			ax[q]   = ax[q+1];
        			ax[q+1] = swap;
      			}
    		}
  		}
		tmp = ax[j/2];
        if (tmp == xr) {
            tmp -= 0.5;
		}
        if (tmp == xl) {
            tmp += 0.5;
		}

        // find the intersection in median
        float alpha = -INFINITY, belta = INFINITY, alpha_m1, alpha_m2, belta_m1, belta_m2, tmp2;
        i = neg;
        while(i > -1) {
			// calculate y and find bigest
			// m1, m2 store the bigest y's two lines' slope
            tmp2 = ((float)data[i].c - data[i].x * tmp) / data[i].y;
            if (tmp2 > alpha) {
                alpha = tmp2;
                alpha_m1 = (float)data[i].x / (-1 * data[i].y);
                alpha_m2 = alpha_m1;
            } else if(tmp2 == alpha) {
                alpha_m2 = (float)data[i].x / (-1 * data[i].y);
			}
            i = data[i].next;
        }
        i = pos;
        while(i > -1) {
			// calculate y and find smallest
			// m1, m2 store the smallest y's two lines' slope
            tmp2 = ((float)data[i].c - data[i].x * tmp) / data[i].y;
            if(tmp2 < belta) {
                belta = tmp2;
                belta_m1 = (float)data[i].x / (-1 * data[i].y);
                belta_m2 = belta_m1;
            } else if(tmp2 == belta) {
                belta_m2 = (float)data[i].x / (-1 * data[i].y);
			}
            i = data[i].next;
        }

        // find the Xm
        if(alpha <= belta) {
            if(-INFINITY != alpha) {
                if(0 >= alpha_m1 * alpha_m2) {
                    printf("%d\n", rounding(alpha));
                    return 0;
                } else if(0 < alpha_m1) {
                    xr = tmp;
				} else {
                    xl = tmp;
				}
            }
        } else {
			float alpha_b, alpha_s, belta_b, belta_s;
			if (alpha_m1 > alpha_m2) {
				alpha_b = alpha_m1;
				alpha_s = alpha_m2;
			} else {
				alpha_b = alpha_m2;
				alpha_s = alpha_m1;
			}

			if (belta_m1 < belta_m2) {
				belta_s = belta_m1;
				belta_b = belta_m2;
			} else {
				belta_s = belta_m2;
				belta_b = belta_m1;
			}
			
            if(alpha_b < belta_s) {
                xl = tmp;
			} else if(alpha_s > belta_b) {
                xr = tmp;
			} else {
                printf("NA\n");
			}
        }
    }
    // when line is 3, process it directly
    if (-INFINITY == xl || INFINITY == xr) {
        printf("-INF\n");
        return 0;
    }
    // concat all line in neg array
    i = neg;
    while (-1 < i && -1 < data[i].next) {
        i = data[i].next;
	}
    data[i].next = pos;
    i = neg;

    // find all intersection
    float inter_x[6], inter_y[6], min_y = INFINITY, tmp2;
    j = -1;
    i = neg;
    while(-1 < data[i].next) {
        ++j;
		// prevent that three lines are parallel
        inter_x[j] = 0;
        inter_y[j] = ((float)data[i].c) / data[i].y;

        if(intersect2(&data[i], &data[data[i].next], &tmp, &tmp2)) {
            ++j;
            inter_x[j] = tmp;
            inter_y[j] = tmp2;
        }
        i = data[i].next;
    }
    ++j;
    inter_x[j] = 0;
    inter_y[j] = ((float)data[i].c) / data[i].y;
    if(intersect2(&data[i], &data[neg], &tmp, &tmp2)) {
        ++j;
        inter_x[j] = tmp;
        inter_y[j] = tmp2;
    }

    // renew the lowest Y_value
    i=0;
    while(i <= j) {
        k = neg;
        int flag = 1;
        while(-1 < k) {
            if(rounding(inter_x[i] * data[k].x + inter_y[i] * data[k].y) > data[k].c) {
                flag = 0;
                break;
            }
            k = data[k].next;
        }
        if (flag) {
            if(inter_y[i] < min_y) {
                min_y = inter_y[i];
            }
        }
        ++i;
    }
    if(INFINITY == min_y) {
        printf("NA\n");
    } else {
        printf("%d\n", rounding(min_y));
    }
    return 0;
}
