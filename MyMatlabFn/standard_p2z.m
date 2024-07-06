%%
function z = standard_p2z(level)

    p_database   = [1000, 900,  850,  700,  600,  500,  400,  300,   250,   200,   150,   100,     70,  50,    30,    20,    10]*100;  %/ in Pa
    p2z_database = [ 111, 988, 1457, 3012, 4206, 5574, 7185, 9164, 10363, 11784, 13608, 16180,  18495, 20576, 23849, 26481, 31055];      %/ in m   (roughly)

    ind = findismember_loop(p_database, level);
    z = p2z_database(ind);

    if length(z) ~= length(level)
        error('The requested level has no corresponding altitude stored in the function!')
    end
end