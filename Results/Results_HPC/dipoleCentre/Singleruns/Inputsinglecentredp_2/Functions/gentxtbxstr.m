function out = gentxtbxstr(Settings)
out = {['fus = ',num2str(Settings{2}*1e-6),' MHz']
    ['mean I = ',num2str(Settings{4}),' µA']
    ['std I = ',num2str(Settings{6}),' µA']
    ['Thermal Noise = ',num2str(Settings{8}),' µV']
    ['dpDistribution = ',Settings{10}]
    ['dpOrientation = ',Settings{12}]
    ['dpI_{space} = ',Settings{14}]};    
    
end