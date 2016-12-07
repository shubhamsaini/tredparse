import React from 'react';
import { Meteor } from 'meteor/meteor';
import { ReactMeteorData } from 'meteor/react-meteor-data';
import { Random } from 'meteor/random';
import Loading from './Loading';
import { ButtonToolbar, Button, FormControl, Form, FormGroup, HelpBlock, Panel } from 'react-bootstrap';
import Documents from '../../api/documents/documents';

const Treds = require('../../api/documents/treds.json');

const Default = {
  httpBAM: '',
  s3BAM: 's3://hli-processed/processed/341087/isaac_align/164648328_S1.bam',
  tred: 'HD',
};

const FormInput = React.createClass({
  propTypes: {
    clickHandler: React.PropTypes.func.isRequired,
  },

  mixins: [ReactMeteorData],

  getMeteorData() {
    let data = {};
    const currentId = this.state.currentId;
    const handle = Meteor.subscribe('documents.view', currentId);
    if (handle.ready()) {
      data.post = Documents.findOne({ _id: currentId });
    }
    return data;
  },

  getInitialState() {
    return {
      currentId: '',
      bam: Default.s3BAM,
      tred: Default.tred,
    };
  },

  handleChange(e) {
    this.setState({ value: e.target.value });
  },

  handleSubmit(e) {
    const currentId = Random.id();
    Meteor.call('shell', { _id: currentId, cmd: 'sleep 1 && pwd' }, (err) => {
      this.setState({ currentId });
    });
  },

  getContent() {
    return (
      <div>
        Current session Id: { this.state.currentId }
        <br />
        Command: { this.data.post.title }
        <br />
        Stdout: { this.data.post.body }
      </div>
    );
  },

  render() {
    const Buttons = Object.keys(Treds).map((b) => {
      return (
        <Button key={ b } onClick={ this.props.clickHandler.bind(null, b) }>
          { b }
        </Button>
      );
    });

    return (
      <Form>
        <FormGroup
          controlId="formBasicText"
        >
          <Panel header={ <strong>BAM file</strong> }>
            <FormControl
              bsSize="sm"
              type="text"
              ref='bam'
              value={ this.state.bam }
              placeholder="Enter sample BAM here"
              onChange={ this.handleChange }
            />
            <FormControl.Feedback />
            <HelpBlock>
                BAM file could be either on <Button bsSize='small' bsStyle='link' onClick={ () => {
                  this.setState({ bam: Default.httpBAM });
                }}>
                  HTTP</Button> or <Button bsSize='small'
                bsStyle='link' onClick={ () => {
                  this.setState({ bam: Default.s3BAM });
                }}>S3</Button>
            </HelpBlock>
          </Panel>
        </FormGroup>

        <FormGroup>
          <Panel header={ <strong>STR locus</strong> }>
            <ButtonToolbar>
              { Buttons }
            </ButtonToolbar>
          </Panel>
        </FormGroup>

        <Button bsStyle='danger' bsSize='large' onClick={ this.handleSubmit }>
          Submit
        </Button>
        <br /><br />
        { this.data.post ? this.getContent() :
          ( this.state.currentId ? <Loading /> : '' )}
      </Form>
    );
  },
});

export default FormInput;
