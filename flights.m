# aux_Airline scheduling problem (aux_ASP)

###############################################################
#                          Constants                          #
###############################################################

# Cities
n = 4;

# Demand
demand = [
     0   450     0     0;
   600     0   450   760;
     0   500     0     0;
     0   700     0     0
];

# Fare
fare = [
     0   175     0     0;
   175     0   230   200;
     0   230     0     0;
     0   200     0     0
];

# Distance
distance = [
     0   260     0     0;
   260     0   375   310;
     0   375     0     0;
     0   310     0     0
];

# Planes
#  cap   spe   pri   uti
planes = [
    50   400  1850    13;
   100   425  3800    12
];

###############################################################
#                       Deduced values                        #
###############################################################

# Air time
time_1 = distance ./ planes(1,2);
time_2 = distance ./ planes(2,2);

time = zeros(size(time_1,1), size(time_1,2), 2);

time(:,:,1) = time_1;
time(:,:,2) = time_2;

###############################################################

# Costs
cost_1 = time(:,:,1) .* planes(1,3);
cost_2 = time(:,:,2) .* planes(2,3);

costs = zeros(size(cost_1,1), size(cost_1,2), 2);

costs(:,:,1) = cost_1;
costs(:,:,2) = cost_2;

###############################################################

# Utilization
utilization = planes(:,4);

###############################################################

# Capacity
capacity = planes(:,1);

###############################################################
#                        User defined                         #
###############################################################

# Minimum number of flights
min_flights = [0 1 0 0; 1 0 1 1; 0 1 0 0; 0 1 0 0];

###############################################################
###############################################################
#                           Solver                            #
###############################################################
###############################################################

# GLPK

flatten = costs(:)';
flatten(~isfinite(flatten)) = 0;
flatten_costs = flatten(flatten > 0);

flatten = demand(:)';
flatten(~isfinite(flatten)) = 0;
flatten_dem = flatten(flatten > 0);

flatten = min_flights(:)';
flatten(~isfinite(flatten)) = 0;
flatten_mins = flatten(flatten > 0);

A = [];
aux_b = vertcat(flatten_dem, flatten_mins)'(:)';
ctype = {};

###############################################################
#                     Demand constraints                      #
###############################################################

for i = 1:size(flatten_dem,2)
    aux = zeros(1, size(flatten_costs,2));
    aux(i) = capacity(1);
    aux(i + size(flatten_dem,2)) = capacity(2);
    A = vertcat(A, aux);
    ctype = vertcat(ctype, 'L');
endfor;

###############################################################
#            Minimum number of flights constraints            #
###############################################################

for i = 1:size(flatten_mins,2)
    aux = zeros(1, size(flatten_costs,2));
    aux(i) = 1;
    aux(i + size(flatten_dem,2)) = 1;
    A = vertcat(A, aux);
    ctype = vertcat(ctype, 'L');
endfor;

###############################################################
#               Aircraft usage time constraints               #
###############################################################

flatten = time(:,:,1)(:)';
flatten(~isfinite(flatten)) = 0;
flatten_time_1 = flatten(flatten > 0);

flatten_time_1 = vertcat(
    flatten_time_1,
    zeros(1, size(flatten_time_1, 2)))'(:)';

flatten = time(:,:,2)(:)';
flatten(~isfinite(flatten)) = 0;
flatten_time_2 = flatten(flatten > 0);

flatten_time_2 = vertcat(
    zeros(1, size(flatten_time_2, 2)),
    flatten_time_2)'(:)';

A = vertcat(A, flatten_time_1);
A = vertcat(A, flatten_time_2);

###############################################################
#                       Loop constants                        #
###############################################################

lower_boundary = zeros(1, size(flatten_costs,2));
upper_boundary = Inf(1, size(flatten_costs,2));
vtype = repmat('I', 1, size(flatten_costs,2));
ctype = vertcat(ctype, 'U');
ctype = vertcat(ctype, 'U');
ctype = cell2mat(ctype);

permutation = [];
min_costs = Inf;
aircraft_x = [];
aircraft_y = [];

###############################################################
#             Find best combinations of aircrafts             #
###############################################################

result = [];
brute_profit = sum(sum(fare .* demand))

for x = 0:3
    for y = 0:3
        needed_utilization = planes(:,4) .* [x y]';

        b = aux_b;
        b(end + 1) = needed_utilization(1);
        b(end + 1) = needed_utilization(2);

        sense = 1;

        [xopt,zmx] = glpk(
            flatten_costs,
            A,
            b,
            lower_boundary,
            upper_boundary,
            ctype,
            vtype,
            sense
        );

        if isfinite(zmx)
            result = vertcat(result,[(brute_profit - zmx) x y]);
            if zmx < min_costs
                min_costs = zmx;
                permutation = xopt;
            endif;
        endif;
    endfor;
endfor;
