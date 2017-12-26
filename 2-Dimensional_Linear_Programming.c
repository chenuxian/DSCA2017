#include <stdio.h>
#include <math.h>

int line_num, pos = -1, neg = -1;

struct Line {
    int a;
    int b;
    int c;
    int prev_i;
    int next_i;
};

int find_intersection(struct Line *l1, struct Line *l2, float *tmp_x) {
    int delta = (l1 -> a * l2 -> b) - (l1 -> b * l2 -> a);
    int delta_x = (l1 -> c * l2 -> b) - (l1 -> b * l2 -> c);
    int delta_y = (l1 -> a * l2 -> c) - (l1 -> c * l2 -> a);

    if (delta != 0) { // unique solution
        *tmp_x = (float)delta_x / delta;
        return 0;
    } else if (delta_x == 0 && delta_y == 0) { // infinite solution i.e. same line
        return 1;
    } else { // no soution i.e. parallel
        return 2;
    }
}

int find_intersection2(struct Line *l1, struct Line *l2, float *tmp_x, float *tmp_y) {
    int delta = (l1 -> a * l2 -> b) - (l1 -> b * l2 -> a);
    int delta_x = (l1 -> c * l2 -> b) - (l1 -> b * l2 -> c);
    int delta_y = (l1 -> a * l2 -> c) - (l1 -> c * l2 -> a);

    if (delta != 0) { // unique solution
        *tmp_x = (float)delta_x / delta;
        *tmp_y = (float)delta_y / delta;
        return 0;
    } else if (delta_x == 0 && delta_y == 0) { // infinite solution i.e. same line
        return 1;
    } else { // no soution i.e. parallel
        return 2;
    }
}

void drop_pos_line(struct Line *line, int i) {
	line_num--;
	if (line[i].next_i != -1 && line[i].prev_i != -1) {
		line[line[i].next_i].prev_i = line[i].prev_i;
		line[line[i].prev_i].next_i = line[i].next_i;
	} else if (line[i].next_i == -1 && line[i].prev_i == -1) {
		pos = -1;
	} else {
		if (line[i].next_i == -1) {
			line[line[i].prev_i].next_i = -1;
		} else {
			pos = line[i].next_i;
			line[line[i].next_i].prev_i = -1;
		}
	}
}

void drop_neg_line(struct Line *line, int i) {
	line_num--;
	if (line[i].next_i != -1 && line[i].prev_i != -1) {
		line[line[i].next_i].prev_i = line[i].prev_i;
		line[line[i].prev_i].next_i = line[i].next_i;
	} else if (line[i].next_i == -1 && line[i].prev_i == -1) {
		neg = -1;
	} else {
		if (line[i].next_i == -1) {
			line[line[i].prev_i].next_i = -1;
		} else {
			neg = line[i].next_i;
			line[line[i].next_i].prev_i = -1;
		}
	}
}

void swap(float *a,float *b) {
	if (*a != *b) {
		*(int *)a ^= *(int *)b;
		*(int *)b ^= *(int *)a;
		*(int *)a ^= *(int *)b;
	}
}

void sort(float a[], int n) {
	int i, j, temp;
	for (i = 0; i < n - 1; ++i) {
    	for(j = 0; j < n - i - 1; ++j) {
      		if(a[j] > a[j + 1]) {
            	swap(&a[j],&a[j + 1]);
			}
      	}
    }
}

int rounding(float a) {
    return (a >= 0) ? (int)(a + 0.5555555555555555) : (int)(a - 0.55555555555555555);
}

int main() {
    int i, count, result;
    int ver_pos = -1, ver_neg = -1, prev1, prev2;
    scanf("%d", &line_num);
    struct Line line[line_num];

    for (i = 0; i < line_num; ++i) {
        scanf("%d", &line[i].a);
        scanf("%d", &line[i].b);
        scanf("%d", &line[i].c);
        line[i].prev_i = -1;
        line[i].next_i = -1;

        // split to slope positive and slope negative
        // if exists vertical line, only leave the critical lines, and then use for boundary
        if (line[i].b == 0) {
			line_num--;
            if (line[i].a > 0) { // x <= c/a
                if (ver_neg == -1) {
                    ver_neg = i;
                } else {
                    if ((float)line[i].c / line[i].a < (float)line[ver_neg].c / line[ver_neg].a) {
                        ver_neg = i;
                    }
                }
            } else { // x >= c/a
                if (ver_pos == -1) {
                    ver_pos = i;
                } else {
                    if ((float)line[i].c / line[i].a > (float)line[ver_pos].c / line[ver_pos].a) {
                        ver_pos = i;
                    }
                }
            }

            if (ver_pos != ver_neg) { // both not -1
                if ((float)line[ver_pos].c / line[ver_pos].a > (float)line[ver_neg].c / line[ver_neg].a) {
                    printf("NA\n");
                    return 0;
                }
            }
        } else {
            if (line[i].b > 0) { // y <= ...
                if (neg == -1) {
                    neg = i;
                } else {
                    line[i].prev_i = prev1;
                    line[prev1].next_i = i;
                }
                prev1 = i;
            } else { // y >= ...
                if (pos == -1) {
                    pos = i;
                } else {
                    line[i].prev_i = prev2;
                    line[prev2].next_i = i;
                }
                prev2 = i;
            }
        }
    } // for (i = 0; i < line_num; ++i)
	
	// all negative lines
	if (pos == -1) {
		printf("-INF\n");
		return 0;
	}

    float xr = INFINITY, xl = -INFINITY, intersection[line_num], tmp_x, tmp_y;

    if (ver_neg != -1) {
        xr = (float)line[ver_neg].c / line[ver_neg].a * -1;
    }
    if (ver_pos != -1) {
        xl = (float)line[ver_pos].c / line[ver_pos].a * -1;
    }

    while (line_num > 3) {
        count = 0;
        // find I+ intersection
        i = pos;
        while (i != -1 && line[i].next_i != -1) {
            result = find_intersection(&line[i], &line[line[i].next_i], &tmp_x);
            if (result == 2) { // parallel
                if (line[i].c > line[line[i].next_i].c) {
					drop_pos_line(line, line[i].next_i);
                } else {
					drop_pos_line(line, i);
                    i = line[i].next_i;
                }
            } else if (result == 1) { // same line
				drop_pos_line(line, line[i].next_i);
            } else {
                if (tmp_x > xr) {
                    // remove bigger slope(-a/b) line
                    if ((double)line[i].a / line[i].b > (double)line[line[i].next_i].a / line[line[i].next_i].b) {
						drop_pos_line(line, line[i].next_i);
                    } else {
						drop_pos_line(line, i);
                        i = line[i].next_i;
                    }
                } else if (tmp_x < xl) {
                    // remove smaller slope(-a/b) line
                    if ((double)line[i].a / line[i].b > (double)line[line[i].next_i].a / line[line[i].next_i].b) {
                        drop_pos_line(line, i);
                        i = line[i].next_i;
                    } else {
						drop_pos_line(line, line[i].next_i);
                    }
                } else {
                    intersection[count++] = tmp_x;
                    i = line[i].next_i;
                }
            }
            i = line[i].next_i;
        }

        // find I- intersection
        i = neg;
        while (i != -1 && line[i].next_i != -1) {
            result = find_intersection(&line[i], &line[line[i].next_i], &tmp_x);
            if (result == 2) { // parallel
                if (line[i].c > line[line[i].next_i].c) {
					drop_neg_line(line, i);
                    i = line[i].next_i;
                } else {
					drop_neg_line(line, line[i].next_i);
                }
            } else if (result == 1) { // same line
				drop_neg_line(line, line[i].next_i);
            } else {
                if (tmp_x > xr) {
                    // remove bigger slope(-a/b) line
                    if ((double)line[i].a / line[i].b < (double)line[line[i].next_i].a / line[line[i].next_i].b) {
						drop_neg_line(line, line[i].next_i);
                    } else {
						drop_neg_line(line, i);
                        i = line[i].next_i;
                    }
                } else if (tmp_x < xl) {
                    // remove smaller slope(-a/b) line
                    if ((double)line[i].a / line[i].b < (double)line[line[i].next_i].a / line[line[i].next_i].b) {
						drop_neg_line(line, i);
                        i = line[i].next_i;
                    } else {
						drop_neg_line(line, line[i].next_i);
                    }
                } else {
                    intersection[count++] = tmp_x;
                    i = line[i].next_i;
                }
            }
            i = line[i].next_i;
        }

        // intersection's median
		sort(intersection, count);
        tmp_x = intersection[(count + 1) / 2 - 1];
		if(tmp_x == xr) {
            tmp_x -= 0.5;
		}
        if(tmp_x == xl) {
            tmp_x += 0.5;
		}
        
		// find the lowest I- and highest I+
		float alpha = -INFINITY, belta = INFINITY, tmp_slope;
		struct Slope {
			float slope;
			int index;
		};
		struct Slope Smax = {-INFINITY, -1}, Smin = {INFINITY, -1}, Tmax = {-INFINITY, -1}, Tmin = {INFINITY, -1};
		i = neg;
		while (i != -1) {
			tmp_y = ((float)line[i].c - line[i].a * tmp_x) / line[i].b;
			if (tmp_y < belta) {
				belta = tmp_y;
				Tmax.slope = Tmin.slope = (float)line[i].a / line[i].b * -1;
				Tmax.index = Tmin.index = i;
			} else if (tmp_y == belta) {
				tmp_slope = (float)line[i].a / line[i].b * -1;
				if (tmp_slope > Tmax.slope) {
					Tmax.slope = tmp_slope;
					Tmax.index = i;
				} else if (tmp_slope < Tmin.slope) {
					Tmin.slope = tmp_slope;
					Tmin.index = i;
				}
			}
			i = line[i].next_i;
		}
		i = pos;
		while (i != -1) {
			tmp_y = ((float)line[i].c - line[i].a * tmp_x) / line[i].b;
			if (tmp_y > alpha) {
				alpha = tmp_y;
				Smax.slope = Smin.slope = (float)line[i].a / line[i].b * -1;
				Smax.index = Smin.index = i;
			} else if (tmp_y == alpha) {
				tmp_slope = (float)line[i].a / line[i].b * -1;
				if (tmp_slope > Smax.slope) {
					Smax.slope = tmp_slope;
					Smax.index = i;
				} else if (tmp_slope < Smin.slope) {
					Smin.slope = tmp_slope;
					Smin.index = i;
				}
			}
			i = line[i].next_i;
		}
		
		// find solution (6 cases)
		if (alpha <= belta) {
			if (alpha != -INFINITY) {
				if (Smin.slope < 0 && Smax.slope < 0) {
                    //drop_pos_line(line, Smax.index);
					xl = tmp_x;
				} else if (Smin.slope > 0 && Smax.slope > 0) {
                    //drop_pos_line(line, Smin.index);
					xr = tmp_x;
				} else {
					// solution is (tmp_x, alpha)
					printf("%d\n", rounding(alpha));
					return 0;
				}
			} else {
				printf("-INF\n");
				return 0;
			}
		} else {
            if (Smax.slope < Tmin.slope) {
                xl = tmp_x;
            } else if (Smin.slope > Tmax.slope) {
                xr = tmp_x;
            } else {
                printf("NA\n");
                return 0;
            }
        }		
    } // while (line_num > 3)

    // when line_num == 3, find solution directly
	if (xl == -INFINITY || xr == INFINITY || pos == -1) {
        printf("-INF\n");
        return 0;
    }
	
	// concatenate all lines, pos will not be -1
    i = pos;
    while(line[i].next_i != -1) {
        i = line[i].next_i;
	}
    line[i].next_i = neg;
	i = pos;

    // find all intersections
    float x[12], y[12], min_y = INFINITY;
	count = 0;
    while (line[i].next_i != -1) {
		// prevent all parallel
        x[count] = 0;
        y[count] = (float)line[i].c / line[i].b;
        
		if (find_intersection2(&line[i], &line[line[i].next_i], &tmp_x, &tmp_y) == 0) {
            ++count;
            x[count] = tmp_x;
            y[count] = tmp_y;
        }
        i = line[i].next_i;
		++count;
    }
	
	// last and first line's intersection
    x[count] = 0;
    y[count] = (float)line[i].c / line[i].b;
    if (find_intersection2(&line[i], &line[pos], &tmp_x, &tmp_y) == 0) {
        ++count;
        x[count] = tmp_x;
        y[count] = tmp_y;
    }

    // because vertical line use for boundary(if exist), need to check, too
    i = pos;
    while (i != -1) {
        if (ver_neg != -1) {
            ++count;
            tmp_x = (float)line[ver_neg].c / line[ver_neg].a; 
            x[count] = tmp_x;
            y[count] = ((float)line[i].c - line[i].a * tmp_x) / line[i].b;
        }
        if (ver_pos != -1) {
            ++count;
            tmp_x = (float)line[ver_pos].c / line[ver_pos].a; 
            x[count] = tmp_x;
            y[count] = ((float)line[i].c - line[i].a * tmp_x) / line[i].b;
        }
        i = line[i].next_i;
    }

    // find min y
    i = 0;
	int j, flag;
    while (i <= count) {
        j = pos;
		flag = 1;
        // check if satisfy all lines(including xl, xr)
        if (x[i] >= xl && x[i] <= xr) {
            while (j != -1) {
                if (rounding(line[j].a * x[i] + line[j].b * y[i]) > line[j].c) {
				    flag = 0;
				    break;
			    }
                j = line[j].next_i;
            }
            if (y[i] < min_y && flag == 1) {
                min_y = y[i];
		    }
        }
       	++i;
    }
    if(min_y == INFINITY) {
        printf("NA\n");
	} else {
        // maybe the solution satisfy three lines but the real solution is -INF
        flag = 1;
        i = pos;
        tmp_y = min_y - 10;
        while (i != -1) {
            if (line[i].a != 0) {
                tmp_x = ((float)line[i].c - line[i].b * tmp_y) / line[i].a;
                // check if satisfy all lines(including xl, xr)
                if (tmp_x >= xl && tmp_x <= xr) {
                    j = pos;
                    while (j != -1) {
                        if (rounding(line[j].a * tmp_x + line[j].b * tmp_y) > line[j].c){
                            flag = 0;
                            break;
                        }
                        j = line[j].next_i;
                    }
                } else {
                    flag = 0;
                    break;
                }
            }
            if (flag == 0) {
                break;
            }
            i = line[i].next_i;
        }
        if (flag == 0) {
            printf("%d\n", rounding(min_y));	
        } else {
            printf("-INF\n");
        }   
	}
    return 0;
}
