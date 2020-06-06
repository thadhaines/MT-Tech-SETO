%clean zero test
a.one = zeros(4,30);
a.two = ones(4);
a.b.one = ones(3);
a.b.z = zeros(50);
a.b.c.two = ones(3);
a.b.c.one = zeros(40);
a.b.c.d.zero = zeros(30);

whos('a')
clearedVars ={};

[a clearedVars] = cleanZeros(a, clearedVars);
whos('a')