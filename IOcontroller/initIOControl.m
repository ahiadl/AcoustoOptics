function [IOdev] = initIOControl()
deviceString = daq.getDevices;
s = daq.createSession('ni');
s.addDigitalChannel(deviceString.ID, 'Port1/Line4', 'OutputOnly');

IOdev.session = s;
IOdev.device = deviceString;

end

