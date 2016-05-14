import rnamotifs2
import math
from fisher import pvalue

def es(data_s, data_e, data_c):
    num_s, num_e, num_c = rnamotifs2.data.dist["s"], rnamotifs2.data.dist["e"], rnamotifs2.data.dist["c"]
    logs = []
    loge = []
    for a in range(0, 4):
        temp_e = []
        temp_s = []
        for val_s, val_e, val_c in zip(data_s[a], data_e[a], data_c[a]): # there are 4 areas
            f1 = -2 * math.log(pvalue(val_s, val_c, num_s-val_s, num_c-val_c).right_tail)
            f2 = -2 * math.log(pvalue(val_e, val_c, num_e-val_e, num_c-val_c).right_tail)
            temp_s.append(f1)
            temp_e.append(f2)
        logs.append(temp_s)
        loge.append(temp_e)
    return logs, loge

def es_apa(data_s, data_e, data_c):
    num_s, num_e, num_c = rnamotifs2.data.dist["s"], rnamotifs2.data.dist["e"], rnamotifs2.data.dist["c"]
    logs = []
    loge = []
    temp_e = []
    temp_s = []
    for val_s, val_e, val_c in zip(data_s, data_e, data_c): # there are 4 areas
        f1 = -2 * math.log(pvalue(val_s, val_c, num_s-val_s, num_c-val_c).right_tail)
        f2 = -2 * math.log(pvalue(val_e, val_c, num_e-val_e, num_c-val_c).right_tail)
        temp_s.append(f1)
        temp_e.append(f2)
    logs.append(temp_s)
    loge.append(temp_e)
    return logs, loge

def test(comps, genome, motif, rcounts, nums): # rcounts
    print "%s.%s.%s: fisher test on real and perm data" % (comps, genome, motif)
    results = {}
    for rt in ["r1", "r2", "r3"]:
        for event_class in ["s", "e"]:
            val_class = rcounts.get("%s.%s" % (rt, event_class), 0)
            val_control = rcounts.get("%s.%s" % (rt, "c"), 0)
            num_class = nums.get("%s.%s" % (rt, event_class), 0)
            num_control = nums.get("%s.%s" % (rt, "c"), 0)
            val = pvalue(val_class, val_control, num_class-val_class, num_control-val_control).right_tail
            results["%s.%s" % (rt, event_class)] = val

            # information gain (g at the end)
            i1 = (num_class/float(num_class+num_control))*math.log(num_class/float(num_class+num_control), 2)
            i2 = (num_control/float(num_class+num_control))*math.log(num_control/float(num_class+num_control), 2)
            i = -(i1+i2)
            c1_num = max(1, float(val_class + val_control)) # dont allow it to be 0
            if (val_class/c1_num)>0 and val_control/c1_num>0:
                c1_num = -( (val_class/c1_num)*math.log(val_class/c1_num, 2) + (val_control/c1_num)*math.log(val_control/c1_num, 2) )
                c1_num = (val_class+val_control)/float(num_class+num_control) * c1_num
            else:
                c1_num = 0
            c2_num = max(1, float( (num_class-val_class) + (num_control-val_control) ) )  # dont allow it to be 0
            if (num_class-val_class)/c2_num>0 and (num_control-val_control)/c2_num>0:
                c2_num = -( ((num_class-val_class)/c2_num)*math.log((num_class-val_class)/c2_num, 2) + ((num_control-val_control)/c2_num)*math.log((num_control-val_control)/c2_num, 2) )
                c2_num = (num_class+num_control-(val_class+val_control))/float(num_class+num_control) * c2_num
            else:
                c2_num = 0
            g = i - (c1_num + c2_num)

            #print rt, event_class, val_class, num_class
            #print rt, "c", val_control, num_control
            #print

            results["%s.%s.g" % (rt, event_class)] = g
            for p in range(0, rnamotifs2.config.perms):
                val_class = rcounts.get("%s.%s.p%s" % (rt, event_class, p), 0)
                val_control = rcounts.get("%s.%s.p%s" % (rt, "c", p), 0)
                num_class = rnamotifs2.perm.ec_dist[p].get(event_class, 0)
                num_control = rnamotifs2.perm.ec_dist[p].get("c", 0)
                val = pvalue(val_class, val_control, num_class-val_class, num_control-val_control).right_tail
                results["%s.%s.p%s" % (rt, event_class, p)] = val

    test_results = {}
    for rt in ["r1", "r2", "r3"]:
        for event_class in ["s", "e"]:
            pval = results["%s.%s" % (rt, event_class)]
            pemp = [results["%s.%s.p%s" % (rt, event_class, p)] for p in range(0, rnamotifs2.config.perms)]
            g = results["%s.%s.g" % (rt, event_class)]
            test_results["%s.%s" % (rt, event_class)] = [pval, pemp, g]
    return test_results

def rtest(comps, genome, motif, rcounts, nums): # rcounts
    print "%s.%s.%s: fisher test on real and perm data" % (comps, genome, motif)
    val_class = rcounts.get("t", 0)
    val_control = rcounts.get("c", 0)
    #num_class = nums.get("t.all", 0) # v_18
    #num_control = nums.get("c.all", 0) # v_18
    num_class = nums.get("t", 0) # v_17
    num_control = nums.get("c", 0) # v_17
    num_all = float(num_class+num_control)
    val_all = float(val_class+val_control)

    fisher = pvalue(val_class, val_control, num_class-val_class, num_control-val_control).right_tail

    # information gain (g at the end)
    i1 = (num_class/num_all)*math.log(num_class/num_all, 2)
    i2 = (num_control/num_all)*math.log(num_control/num_all, 2)
    i = -(i1+i2)
    c1_num = max(1, val_all) # dont allow it to be 0
    if (val_class/c1_num)>0 and val_control/c1_num>0:
        c1_num = -( (val_class/c1_num)*math.log(val_class/c1_num, 2) + (val_control/c1_num)*math.log(val_control/c1_num, 2) )
        c1_num = val_all/num_all * c1_num
    else:
        c1_num = 0
    c2_num = max(1, float( (num_class-val_class) + (num_control-val_control) ) )  # dont allow it to be 0
    if (num_class-val_class)/c2_num>0 and (num_control-val_control)/c2_num>0:
        c2_num = -( ((num_class-val_class)/c2_num)*math.log((num_class-val_class)/c2_num, 2) + ((num_control-val_control)/c2_num)*math.log((num_control-val_control)/c2_num, 2) )
        c2_num = (num_all-val_all)/num_all * c2_num
    else:
        c2_num = 0
    ig = i - (c1_num + c2_num)

    p_emp = []
    for p in range(0, rnamotifs2.config.perms):
        val_class = rcounts.get("%s.p%s" % (event_class, p), 0)
        val_control = rcounts.get("c.p%s" % p, 0)
        num_class = rnamotifs2.perm.ec_dist[p].get(event_class, 0)
        num_control = rnamotifs2.perm.ec_dist[p].get("c", 0)
        val = pvalue(val_class, val_control, num_class-val_class, num_control-val_control).right_tail
        p_emp.append(val)

    return (fisher, p_emp, ig)
