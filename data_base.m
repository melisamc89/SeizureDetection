
%%estructure patients will contain information from the data base
patients{1}.name='WL';
patients{2}.name='RF';
patients{3}.name='MA';
patients{4}.name='CM';
patients{5}.name='CG';


%% patien1=wagner;
patients{1}.region.name={'InA','InP','Ami','HiC','Hes'};
patients{1}.region.init{2}={[1,2,3],[5,6,7],[],[],[]};
patients{1}.region.init{3}={[4,5,6],[],[],[],[]};
patients{1}.region.init{4}={[],[5,6],[],[],[]};
patients{1}.region.init{5}={[],[6],[],[],[]};
patients{1}.region.init{6}={[3,5,6,7],[],[],[],[]};
patients{1}.region.ienit{7}={[2,3],[5,6],[],[],[]};
patients{1}.region.init{9}={[],[],[],[],[]};
patients{1}.region.init{11}={[],[],[],[],[]};
patients{1}.region.init{12}={[],[],[],[],[]};
patients{1}.region.init{14}={[],[],[],[],[]};


patients{1}.CSS=[0,0,1,5,4,6,1,0,0,0,1,0,0,0];
patients{1}.Cycle=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]; %% 0 for awake, 1 for sleep
patients{1}.side=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]; %% 0 for right, 1 for left
patients{1}.propagation=[0,0,1,0,1,1,1,0,0,0,0,0,0,0]; % 0 for no, 1 for yes

patients{1}.crisis.electrodo_init=[0,1464,208429,40530,380727,804839,14606,0,0,0,0,0,0,0,0];
patients{1}.crisis.clinic_init=[0,1464,255707,157313,423610,905216,28729,0,0,0,0,0,0,0];
patients{1}.crisis.AOC=[0,0,257745,157313,433975,849836,978631,0,0,0,0,0,0,0,0];
patients{1}.crisis.propagation=[0,0,231156,0,458504,849836,22610,0,0,0,0,0,0,0];
patients{1}.crisis.electrodo_end=[0,89250,397738,309182,667242,1183216,137736,0,0,0,0,0,0,0];
patients{1}.crisis.clinic_end=[0,0,0,0,0,0,0,0,0,0,0,0,0,0];

%% patien2=rizzio;

patients{2}.region.name={'CiA','CiP','GRS','HCu','HCz'};
patients{2}.region.init{1}={[],[],[],[],[]};
patients{2}.region.init{2}={[],[7,8,9],[],[],[]};
patients{2}.region.init{3}={[],[7,8,9],[],[],[]};
patients{2}.region.init{4}={[],[],[],[],[]};
patients{2}.region.init{6}={[],[],[],[],[]};


patients{2}.CSS=[6,9,9,7,0,0];
patients{2}.Cycle=[0,0,0,0,0,0]; %% 0 for awake, 1 for sleep
patients{2}.side=[0,1,1,0,0,0]; %% 0 for right, 1 for left
patients{2}.propagation=[0,1,1,0,0,0];

patients{2}.crisis.electrodo_init=[0,821320,626758,0,0,0];
patients{2}.crisis.clinic_init=[0,843743,735277,0,0,0];
patients{2}.crisis.AOC=[0,912806,735277,0,0,0];
patients{2}.crisis.propagation=[0,843743,641141,0,0,0];
patients{2}.crisis.electrodo_end=[0,1031647,880443,0,0,0];
patients{2}.crisis.clinic_end=[0,0,0,0,0,0];

%% patient3 = molina

patients{3}.region.name={'AmD','HkD','HqD','AmI','HkI','HqI'};
patients{3}.region.init{1}={[],[10,11],[19,20],[],[],[]};
patients{3}.region.init{2}={[],[],[],[28,29],[37,38],[]};
patients{3}.region.init{3}={[],[],[],[28,29],[37,38],[46,47]};
patients{3}.region.init{4}={[],[],[10,11],[19,20],[],[]};


patients{3}.CSS=[9,9,7,9];
patients{3}.Cycle=[0,0,0,1]; %% 0 for awake, 1 for sleep
patients{3}.side=[0,1,1,0]; %% 0 for right, 1 for left
patients{3}.propagation=[1,1,1,1];

patients{3}.crisis.electrodo_init=[256129,0,719999,2532782];
patients{3}.crisis.clinic_init=[260202,0,728105,2567873];
patients{3}.crisis.AOC=[260182,0,0,2567873];
patients{3}.crisis.propagation=[272470,0,746580,2545302];
patients{3}.crisis.electrodo_end=[410728,0,923449,2727414];
patients{3}.crisis.clinic_end=[0,0,0,0];  %

%% patient4 = Capelletti

patients{4}.region.name={'LPr','Hca','LPo','GiPo','Hes','GiPr'};
patients{4}.region.init{37}={[],[],[],[],[],[]};
patients{4}.region.init{38}={[],[],[],[],[],[]};
patients{4}.region.init{48}={[],[],[],[],[],[]};
patients{4}.region.init{50}={[],[],[],[],[],[]};
patients{4}.region.init{58}={[],[],[],[],[],[]};
patients{4}.region.init{62}={[],[],[],[],[],[]};
patients{4}.region.init{66}={[],[],[],[],[],[]};
patients{4}.region.init{68}={[],[],[],[],[],[]};
patients{4}.region.init{69}={[],[],[],[],[],[]};

patients{4}.CSS=zeros(1,69);
patients{4}.Cycle=zeros(1,69);
patients{4}.side=zeros(1,69);
patients{4}.propagation=zeros(1,69);

patients{4}.crisis.electrodo_init=zeros(1,69);
patients{4}.crisis.clinic_init=zeros(1,69);
patients{4}.crisis.AOC=zeros(1,69);
patients{4}.crisis.propagation=zeros(1,69);
patients{4}.crisis.electrodo_end=zeros(1,69);
patients{4}.crisis.clinic_end=zeros(1,69);

%% patient5 = Cardona

patients{5}.region.name={'AmI','HkI','HqI','AmD','HkD','HqD'};
patients{5}.region.init{3}={[],[],[],[],[],[]};
patients{5}.region.init{6}={[],[],[],[],[],[]};
patients{5}.region.init{8}={[],[],[],[],[],[]};

patients{5}.CSS=zeros(1,8);
patients{5}.Cycle=zeros(1,8);
patients{5}.side=zeros(1,8);
patients{5}.propagation=zeros(1,8);

patients{5}.crisis.electrodo_init=zeros(1,8);
patients{5}.crisis.clinic_init=zeros(1,8);
patients{5}.crisis.AOC=zeros(1,8);
patients{5}.crisis.propagation=zeros(1,8);
patients{5}.crisis.electrodo_end=zeros(1,8);
patients{5}.crisis.clinic_end=zeros(1,8);

%% save data_base
save('data_patients','patients','-mat')