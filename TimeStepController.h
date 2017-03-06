#pragma once
#define TIME_STEP_CONTROLLER
#define _CRT_SECURE_NO_WARNINGS

#define TIME_ACCURACY 1e-4

class TimeStepController {
public:
	TimeStepController() {}
	TimeStepController(float totalTime, float fps, float defaultDt) {
		_totalTime = totalTime;
		_spf = 1.0f / fps;
		_frameFlag = false;
		_currentTime = 0;
		_frameTime = _spf;
		_frameCnt = 0;
		_defaultDt = defaultDt;
	}
	TimeStepController(int totalFrame, float fps, float defaultDt) {
		_spf = 1.0f / fps;
		_totalTime = totalFrame * _spf;
		_frameFlag = false;
		_currentTime = 0;
		_frameTime = _spf;
		_frameCnt = 0;
		_defaultDt = defaultDt;
	}
	bool isFinished() {
		return _currentTime>_totalTime;
	}
	float getStepDt() {
		float dt = _defaultDt;
		float supposeNextTime = _currentTime + dt;
		if (supposeNextTime - _frameTime>-TIME_ACCURACY) {
			dt = _frameTime - _currentTime;
			_currentTime = _frameTime;
			_frameTime += _spf;
			_frameFlag = true;
		}
		else {
			_currentTime = supposeNextTime;
			_frameFlag = false;
		}
		return dt;
	}
	/* Whenever the dt * fps -> one frame duration,
	the frameFlag is set, indicating a new actual 'frame' shoule be generated.
	 */
	bool isFrameTime(int &frameCnt) {
		if (_frameFlag) {
			frameCnt = _frameCnt++;
		}
		return _frameFlag;
	}

	float _currentTime;

protected:
	float _totalTime;
	float _spf;
	float _defaultDt;
	float _frameTime;
	bool  _frameFlag;
	int _frameCnt;
};