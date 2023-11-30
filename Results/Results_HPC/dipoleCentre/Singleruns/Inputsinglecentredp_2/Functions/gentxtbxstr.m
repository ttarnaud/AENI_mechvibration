function out = gentxtbxstr(Settings)
out = {['fus = ',num2str(Settings{2}*1e-6),' MHz']
    ['mean I = ',num2str(Settings{4}),' �A']
    ['std I = ',num2str(Settings{6}),' �A']
    ['Thermal Noise = ',num2str(Settings{8}),' �V']
    ['dpDistribution = ',Settings{10}]
    ['dpOrientation = ',Settings{12}]
    ['dpI_{space} = ',Settings{14}]};    
    
end