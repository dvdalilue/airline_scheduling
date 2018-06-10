# Airline scheduling problem (ASP)

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
time_1 = planes(1,2) ./ distance;
time_2 = planes(2,2) ./ distance;

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

# Number of aircraft of type k
aircraft = [3,3];

###############################################################
#                           Solver                            #
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
b = vertcat(flatten_dem, flatten_mins)'(:)';
ctype = {};

for i = 1:size(flatten_dem,2)
    aux = zeros(1, size(flatten_costs,2));
    aux(i) = capacity(1);
    aux(i + size(flatten_dem,2)) = capacity(2);
    A = vertcat(A, aux);
    ctype = vertcat(ctype, 'L');
endfor;

# The next should be commented or deleted, it is just for testing
% A = [];

for i = 1:size(flatten_mins,2)
    aux = zeros(1, size(flatten_costs,2));
    aux(i) = 1;
    aux(i + size(flatten_dem,2)) = 1;
    A = vertcat(A, aux);
    ctype = vertcat(ctype, 'L');
endfor;

# The next should be commented or deleted, it is just for testing
% A = [];

needed_utilization = planes(:,4) .* aircraft';

flatten = time(:,:,1)(:)';
flatten(~isfinite(flatten)) = 0;
flatten_time = flatten(flatten > 0);

flatten_time = vertcat(flatten_time, zeros(1, size(flatten_time, 2)))'(:)';

A = vertcat(A, flatten_time);
b(end + 1) = needed_utilization(1);
ctype = vertcat(ctype, 'U');

flatten = time(:,:,2)(:)';
flatten(~isfinite(flatten)) = 0;
flatten_time = flatten(flatten > 0);

flatten_time = vertcat(zeros(1, size(flatten_time, 2)), flatten_time)'(:)';

A = vertcat(A, flatten_time);
b(end + 1) = needed_utilization(2);
ctype = vertcat(ctype, 'U');

sense = 1;
ctype = cell2mat(ctype);

lower_boundary = zeros(1, size(flatten_costs,2));
upper_boundary = Inf(1, size(flatten_costs,2));
vtype = repmat('I', 1, size(flatten_costs,2));

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