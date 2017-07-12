int: sequence_length = 10;

array[0..sequence_length] of var int: t; 
% don't write something like array[1, 1..sequence_length]
% its dimension is one, and MiniZinc knows that

constraint forall(i in 0..sequence_length)(
    t[i] = sum(k in 0..sequence_length)((bool2int(t[k] == i))));

solve satisfy;
% funny thing, if we initiate as array[1..sequence_length]
% then we will always get something like [0,0,0,0,0]
