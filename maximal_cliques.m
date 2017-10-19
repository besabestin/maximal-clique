%Author Bisrat Tekle Aweno in DEJ Technology GMBH
%Rostock, Germany
%License would be free for anyone to copy and change and distribute or use in commercial 
%products given the above copyright information is provided.

%matlab translation of the code found in 
%http://www.dcs.gla.ac.uk/~pat/jchoco/clique/indSetMachrahanish/papers/tomita2006.pdf

%receives the graph in adjacency matrix
%currently using the MCQ function, next improvement will bring the main function to work.

function tomita(M, use_tomita)
	%[R No] = number_sort(R, []);

	%first test case
	% R = {};
	
	% R{1,1} = 1;
	% R{1,2} = [2, 3];
	% R{2,1} = 2;
	% R{2,2} = [1, 3, 4];
	% R{3,1} = 3;
	% R{3,2} = [1, 2, 4, 5];
	% R{4,1} = 4;
	% R{4,2} = [2, 3];
	% R{5,1} = 5;
	% R{5,2} = [3];

	%second test case
	% R = {};
	
	% R{1,1} = 1;
	% R{1,2} = [2, 3];
	% R{2,1} = 2;
	% R{2,2} = [1, 3, 4];
	% R{3,1} = 3;
	% R{3,2} = [1, 2, 4, 6, 7];
	% R{4,1} = 4;
	% R{4,2} = [2, 3, 5, 6, 7];
	% R{5,1} = 5;
	% R{5,2} = [4];
	% R{6,1} = 6;
	% R{6,2} = [3, 4, 7];
	% R{7,1} = 7;
	% R{7,2} = [3, 4, 6];

	%[R No] = number_sort(R, []);
	%Qm

	clear('global');

	R = to_adjacency_list(M);

	if use_tomita
		tomita_maximal_cliques(R);
	else
		MCQ(R);
	end


endfunction

%changes an adjacency matrix
%to an adjacency list
function L = to_adjacency_list(M)
	L = {};
	len = size(M, 1);
	for i = 1:len
		L{i, 1} = i;
		L{i, 2} = [];
		for j = 1:len
			if M(i, j)
				L{i, 2} = [L{i, 2} j];
			end
		end
	end
endfunction

%this is the branch and bound entry for the np-hard problem.
function MCQ(G)
	global Q = {};
	global Q_max = {};

	%disp(['amazing: ', num2str(size(Q_max, 1))]);

	%Q_max

	%------------------------------%
	%%--%%--SORT PART BEGINS--%%--%%
	%------------------------------%

	%sort vertices
	R = {};
	rcount = 1;

	%calculate the degree for all of them.
	for i = 1:size(G, 1)
		deg(i) = length(G{i, 2});
	end

	while(size(R, 1) < size(G, 1))
		[M I] = max(deg);
		if(M == -1)
			break;
		end
		R(rcount, :) = G(I, :);
		rcount =  rcount + 1;
		deg(I) = -1;
	end

	%check if they are in descending order of degrees.


	G = R;

	%so the above part works fine.

	g_deg = get_max_degree(G);
	%disp(['max degree', num2str(g_deg)]);
	fflush(stdout);
	for i = 1:g_deg
		No(G{i, 1}) = i;
	end
	for i = (g_deg + 1):size(G, 1)
		No(G{i, 1}) = g_deg + 1;
	end
	%disp(['the No thing ', num2str(No)]);
	expand(G, No);
endfunction


%TODO: Proof for two kinds of vertices.
function Q_max = tomita_maximal_cliques(G)

	global Q = {};
	global Q_max = {};

	%------------------------------%
	%%--%%--SORT PART BEGINS--%%--%%
	%------------------------------%

	g_deg = get_max_degree(G);

	i = size(G, 1);
	R = G;
	G = {};
	R_min = min_deg_set(R);

	%disp(['the maximum degree is: ', num2str(i), ' ', num2str(g_deg)]);

	disp(['size of R_min', num2str(size(R_min, 1))]);

	while size(R_min, 1) ~= size(R, 1)
		if size(R_min, 1) >= 2
			p = find_min_ex_deg(R_min);
		else
			p = R_min{1,1}; %TODO: what if R_min{1,1} is non-existent
		end

		%disp(['node number: ', num2str(findNode(R_min, p))]);

		V(i, :) = R_min(findNode(R_min, p), :);
		%disp(['remove node: ', num2str(p)]);
		R = remove_node(R, p);
		i = i - 1;
		for j = 1:size(R, 1)
			adj = R{j, 2};
			%use find buddy
			if any(adj == p)
				%TODO: weird thing here
				%the algorithm just does deg = deg - 1
				%for which I just removed the adjacent thing.
				R{j, 2} = R{j, 2}(R{j, 2}~=p);
			end
		end
		R_min = min_deg_set(R);
	end

	%%--%%--SORT PART ENDS--%%--%%

	[R_min, No] = number_sort(R_min, []);
	for j = 1:size(R_min, 1)
		V(j, :) = R_min(j, :);
	end

	%so I is the index of p
	[m I] = max(No);
% 	if m == -1 || size(R_min, 1) == 0
% 		disp('shouldnt be here nononono!');
% 	end
	%find the index in R that corresponds to I
	%q = findNode(R_min, I);

	mmax = size(R_min, 1) + g_deg - m;
	m = m + 1;
	i = size(R_min, 1) + 1;

	gotoActive = false;

	while i <= mmax
		if i > size(V, 1)
			%TODO: THE GOTO THINGIE, for God's sake
			gotoActive = true;
			break;
		else
			No(V{i, 1}) = m;
			m = m + 1;
			i = i + 1;
		end
	end

	if !gotoActive
		for i = (mmax+1):size(V, 1)
			%disp(['sayin what? ', num2str(i)]);
			No(V{i, 1}) = g_deg + 1;
		end
	end

	%TODO: NOT SURE IF IT WAS ACTUALLY HERE.
	respects = true;
	for i = 1:size(R_min, 1)
		if length(R_min{i, 2}) ~= (size(R_min, 1) - 1)
			respects = false;
			break;
		end
	end

	if respects
		Q_max = R_min;
	end

	expand(V, No);
endfunction


function maxdeg = get_max_degree(G)
	maxdeg = -1;
	for i = 1:size(G, 1)
		%disp(['length: ', num2str(length(G{i, 2}))]);
		maxdeg = max(length(G{i, 2}), maxdeg);
	end
endfunction

%set of vertices with the minimum degree in R
function R_min = min_deg_set(R)
	minI = -1;
	minDeg = inf;
	for i = 1:size(R, 1)
		if length(R{i, 2}) < minDeg
			minI = i;
			minDeg = length(R{i, 2});
		elseif length(R{i, 2}) == minDeg
			minI = [minI i];
		end
	end
	R_min = {};
	%doesn't make sense to maintain the name of the vertex
	R_min(1:length(minI),:) = R(minI,:);
endfunction

function I = find_min_ex_deg(R)
	%find minimum ex degree
	mindeg = inf;
	for i = 1:size(R, 1)
		thisdeg = ex_deg(R, R{i, 1});
		if(mindeg > thisdeg)
			mindeg = thisdeg;
			I = R{i, 1};
		end
	end
endfunction

function deg = ex_deg(R, p)
	deg = 0;
	I = findNode(R, p);
	if I == -1
		return;
	end
	adj = R{I, 2};
	deg = length(adj);
	for i = 1:length(adj)
		j = findNode(R, adj(i));
		if j ~= -1
			adjInner = R{j, 2};
			deg = deg + length(adjInner);
		end
	end
endfunction

%returns the index of the node in the node array
function Ip = findNode(R, val)
	%disp(['the val', num2str(val)]);
	Ip = -1;
	for j = 1:size(R, 1)
		%disp(['talk', num2str(R{j, 1}), ' ', num2str(val)]);
		if R{j,1} == val
			%disp('here');
			Ip = j;
			break;
		end
	end
endfunction

function P = remove_node(Q, val)
	i = findNode(Q, val);
	P = {};
	if i ~= 1
		P(1:(i - 1), :) = Q(1:(i - 1), :);
	end
	qlen = size(Q, 1);
	if qlen ~= 0
		P(i:(qlen - 1), :) = Q((i+1):qlen, :);
	end
endfunction

function [R No] = expand(R, No)
	%then probably this specific No value has to be set to negative
	%Q = {}; %I think these are global variables
	%Q_max = {};
	global Q;
	global Q_max;
	while size(R, 1) ~= 0
		%so I is the index of p
		[d I] = max(No);
		if d == -1 || size(R, 1) == 0
			break;
		end
		%find the index in R that corresponds to I
		Ip = findNode(R, I);
		if Ip == -1
			break;
		end

		p = R(Ip, :);

		Rp = {};

		if (size(Q, 1) + No(p{1,1})) > size(Q_max, 1)
			%union operation
			%shouldn't repeat elements
			unI = findNode(Q, p{1,1});
			if unI == -1
				Q(size(Q, 1) + 1, :) = p; 
				%disp(['added to Q: ', num2str(p{1,1}), ': ', num2str(p{1,2})]);
				%Q
				%otherwise it is already there.Q
			end
			adj = p{1,2};
			cnt = 1;
			%disp('testing');
			%R
			%disp(['against: ', num2str(adj)]);
			for k = 1:length(adj)
				adjI = findNode(R, adj(k));
				if adjI ~= -1
					Rp(cnt,:) = R(adjI, :); %todo: but only the ones that are actually there.
					cnt = cnt + 1;
				end
			end
			%disp(['count: ', num2str(cnt), ' q size: ', num2str(size(Q, 1)), ' qmax size: ', num2str(size(Q_max, 1))]);
			if cnt > 1
				[Rp No1] = number_sort(Rp, []);
				[Rp No1] = expand(Rp, No1);
			elseif size(Q, 1) > size(Q_max, 1)
				Q_max = Q;
				%Q_max
			end
			%remove p from Q
			Q = remove_node(Q, p{1,1});
		else
			return;
		end

		No(p{1,1}) = -1;
		%remove the node
		R = remove_node(R, p{1,1});
	end
endfunction


%function seems clear now.
function [R No] = number_sort(R, No)
	%R = {};
	
	% R{1,1} = 1;
	% R{1,2} = [2,3,4];
	% R{2,1} = 2;
	% R{2,2} = [1,3,4];
	% R{3,1} = 3;
	% R{3,2} = [1,2,4];

	%the coloring part
	maxno = 0;
	C = {};
	rIndex = 1;
	cCount = [];

	while rIndex <= size(R, 1)
		
		p = R(rIndex, :); %R is ...
		k = 1;

		while tomita_intersection(C, k, p)
			k = k + 1;
		end

		if k > maxno
			maxno = k;
			%C{maxno, 1} =  {};
		end
		
		No(p{1,1}) = k;

		if(length(cCount) < k)
			%len = 1;
			cCount(k) = 1;
		else
			cCount(k) = cCount(k) + 1;
		end

		C{k, cCount(k)} = p;

		%forget R = R - {p}
		rIndex = rIndex + 1;

	end

	%C

	R = {};%redo R here.

	%the sorting part
	i = 1;
	for k = 1:maxno
		for j = 1:size(C(k,:), 2)
			%disp(['size of (C(k,:), 2) ', num2str(size(C(k,:), 2))])
			if(size (C{k, j}, 1) ~= 0)
				R(i, :) = C{k, j};
				i = i + 1;
			end
		end
	end

	%disp(['the No array: ', num2str(No)]);
	%the No array is fine as well.


endfunction

function theyIntersect = tomita_intersection(C, k, p)
	if size(C, 1) < k
		theyIntersect = false;
		return;
	end
	for i = 1:length(p{1,2})
		for j = 1:size(C(k, :), 2)
			if size(C{k, j}, 1) ~= 0 && p{1,2}(i) == C{k,j}{1,1}
				theyIntersect = true;
				return;
			end
		endfor
	endfor
	theyIntersect = false;
endfunction
