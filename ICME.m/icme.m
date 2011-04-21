
clear all;
configPath = '../data/config';
disp(['Reading config from file: ', configPath]);
config = Config(configPath);
quality = 'g';
disp(['Filtering events by quality: ', quality])
config.read().filter(quality);
disp(['Number of events to analyze: ', num2str(config.count)]);
disp('The list of events:');
for iEvent = 1:config.count
    disp(['    ', datestr(config.data(iEvent).beginDate, 'yyyy-mm-dd HH:MM')]);
end
for iEvent = 1:config.count
    disp(['Analyzing event ', datestr(config.data(iEvent).beginDate, ...
        'yyyy-mm-dd HH:MM')]);
    dataPath = ['../data/', Config.getDataFilename(...
        config.data(iEvent).spacecraft, config.data(iEvent).samplingInterval)];
    disp(['Reading data from file: ', dataPath]);

    beginDate = addtodate(config.data(iEvent).beginDate, -1, 'day');
    endDate = addtodate(config.data(iEvent).endDate, 1, 'day');

    data = Data(dataPath);
    data.read(beginDate, endDate);
    event = Event();
    event.data = data.copy().filter(config.data(iEvent).beginDate, ...
        config.data(iEvent).endDate);
    event.config = config.data(iEvent);

    if event.config.toComputeGsr
        disp('Performing GSR analysis...');
        event.computeGsr();
        disp('De Hoffmann Teller frame speed:');
        disp(['V = ', num2str(event.gsr.dht.Vx), '; ', ...
            num2str(event.gsr.dht.Vy), '; ', num2str(event.gsr.dht.Vz)]);
        disp(['correlation coefficient = ', num2str(event.gsr.dht.cc)]);
    end

    if event.config.toComputeMva
        disp('Performing MVA analysis...');
        event.computeMva();
        disp('Results of MVAB:');
        disp(['x = ', num2str(event.mva.mvab.x(1)), '; ', ...
            num2str(event.mva.mvab.x(2)), '; ', num2str(event.mva.mvab.x(3))]);
        disp(['y = ', num2str(event.mva.mvab.y(1)), '; ', ...
            num2str(event.mva.mvab.y(2)), '; ', num2str(event.mva.mvab.y(3))]);
        disp(['z = ', num2str(event.mva.mvab.z(1)), '; ', ...
            num2str(event.mva.mvab.z(2)), '; ', num2str(event.mva.mvab.z(3))]);
        disp(['criterion = ', num2str(event.mva.mvab.criterion)]);
        disp('Results of MVUB:');
        disp(['x = ', num2str(event.mva.mvub.x(1)), '; ', ...
            num2str(event.mva.mvub.x(2)), '; ', num2str(event.mva.mvub.x(3))]);
        disp(['y = ', num2str(event.mva.mvub.y(1)), '; ', ...
            num2str(event.mva.mvub.y(2)), '; ', num2str(event.mva.mvub.y(3))]);
        disp(['z = ', num2str(event.mva.mvub.z(1)), '; ', ...
            num2str(event.mva.mvub.z(2)), '; ', num2str(event.mva.mvub.z(3))]);
        disp(['criterion = ', num2str(event.mva.mvub.criterion)]);
    end
end

