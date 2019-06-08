deviceString = daq.getDevices;
s = daq.createSession('ni');
s.addDigitalChannel(deviceString.ID, 'Port1/Line4', 'inputOnly');
while(true)
    data = inputSingleScan(s)
end