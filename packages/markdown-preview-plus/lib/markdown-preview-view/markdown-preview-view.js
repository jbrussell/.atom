"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
const path = require("path");
const atom_1 = require("atom");
const lodash_1 = require("lodash");
const fs = require("fs");
const renderer = require("../renderer");
const markdownIt = require("../markdown-it-helper");
const util_1 = require("../util");
const util = require("./util");
const webview_handler_1 = require("./webview-handler");
const image_watch_helper_1 = require("../image-watch-helper");
const pdf_export_util_1 = require("./pdf-export-util");
const macros_util_1 = require("../macros-util");
class MarkdownPreviewView {
    constructor(renderLaTeX = util_1.atomConfig().mathConfig
        .enableLatexRenderingByDefault) {
        this.renderLaTeX = renderLaTeX;
        this.emitter = new atom_1.Emitter();
        this.disposables = new atom_1.CompositeDisposable();
        this.destroyed = false;
        this.loading = true;
        this.changeHandler = () => {
            util_1.handlePromise(this.renderMarkdown());
            const pane = atom.workspace.paneForItem(this);
            if (pane !== undefined && pane !== atom.workspace.getActivePane()) {
                pane.activateItem(this);
            }
        };
        this.initialRenderPromise = new Promise((resolve) => {
            this.handler = new webview_handler_1.WebviewHandler(async () => {
                const config = util_1.atomConfig();
                await this.handler.init({
                    userMacros: macros_util_1.loadUserMacros(),
                    mathJaxConfig: config.mathConfig,
                    context: 'live-preview',
                });
                await this.handler.setBasePath(this.getPath());
                this.emitter.emit('did-change-title');
                resolve(this.renderMarkdown().then(() => {
                    this.loading = false;
                }));
            });
            this.runJS = this.handler.runJS.bind(this.handler);
            this.imageWatcher = new image_watch_helper_1.ImageWatcher(this.handler.updateImages.bind(this.handler));
            MarkdownPreviewView.elementMap.set(this.element, this);
        });
        this.handleEvents();
    }
    get element() {
        return this.handler.element;
    }
    static viewForElement(element) {
        return MarkdownPreviewView.elementMap.get(element);
    }
    destroy() {
        if (this.destroyed)
            return;
        this.destroyed = true;
        this.imageWatcher.dispose();
        this.disposables.dispose();
        this.handler.destroy();
        MarkdownPreviewView.elementMap.delete(this.element);
    }
    onDidChangeTitle(callback) {
        return this.emitter.on('did-change-title', callback);
    }
    onDidChangeMarkdown(callback) {
        return this.emitter.on('did-change-markdown', callback);
    }
    toggleRenderLatex() {
        this.renderLaTeX = !this.renderLaTeX;
        this.changeHandler();
    }
    getDefaultLocation() {
        return util_1.atomConfig().previewConfig.previewDock;
    }
    getIconName() {
        return 'markdown';
    }
    getSaveDialogOptions() {
        let defaultPath = this.getPath();
        if (defaultPath === undefined) {
            const projectPath = atom.project.getPaths()[0];
            defaultPath = 'untitled.md';
            if (projectPath) {
                defaultPath = path.join(projectPath, defaultPath);
            }
        }
        defaultPath += '.' + util_1.atomConfig().saveConfig.defaultSaveFormat;
        return { defaultPath };
    }
    saveAs(filePath) {
        if (filePath === undefined)
            return;
        if (this.loading)
            throw new Error('Preview is still loading');
        const { name, ext } = path.parse(filePath);
        if (ext === '.pdf') {
            util_1.handlePromise(this.getMarkdownSource().then(async (mdSource) => pdf_export_util_1.saveAsPDF(mdSource, this.getPath(), this.getGrammar(), this.renderLaTeX, filePath)));
        }
        else {
            util_1.handlePromise(this.getHTMLToSave(filePath).then(async (html) => {
                const fullHtml = util.mkHtml(name, html, this.renderLaTeX, await this.handler.getTeXConfig());
                fs.writeFileSync(filePath, fullHtml);
                return atom.workspace.open(filePath);
            }));
        }
    }
    didScrollPreview(_min, _max) {
    }
    openSource(initialLine) {
        const path = this.getPath();
        if (path === undefined)
            return;
        util_1.handlePromise(atom.workspace.open(path, {
            initialLine,
            searchAllPanes: true,
        }));
    }
    syncPreview(line, flash) {
        util_1.handlePromise(this.handler.sync(line, flash));
    }
    openNewWindow() {
        const path = this.getPath();
        if (!path) {
            atom.notifications.addWarning('Can not open this preview in new window: no file path');
            return;
        }
        atom.open({
            pathsToOpen: [`markdown-preview-plus://file/${path}`],
            newWindow: true,
        });
        util.destroy(this);
    }
    handleEvents() {
        this.disposables.add(atom.grammars.onDidAddGrammar(() => lodash_1.debounce(() => {
            util_1.handlePromise(this.renderMarkdown());
        }, 250)), atom.grammars.onDidUpdateGrammar(lodash_1.debounce(() => {
            util_1.handlePromise(this.renderMarkdown());
        }, 250)), atom.commands.add(this.element, {
            'core:move-up': () => this.element.scrollBy({ top: -10 }),
            'core:move-down': () => this.element.scrollBy({ top: 10 }),
            'core:copy': () => {
                util_1.handlePromise(this.copyToClipboard());
            },
            'markdown-preview-plus:open-dev-tools': () => {
                this.handler.openDevTools();
            },
            'markdown-preview-plus:new-window': () => {
                this.openNewWindow();
            },
            'markdown-preview-plus:print': () => {
                this.handler.print();
            },
            'markdown-preview-plus:zoom-in': () => {
                this.handler.zoomIn();
            },
            'markdown-preview-plus:zoom-out': () => {
                this.handler.zoomOut();
            },
            'markdown-preview-plus:reset-zoom': () => {
                this.handler.resetZoom();
            },
            'markdown-preview-plus:sync-source': async (_event) => {
                const line = await this.handler.syncSource();
                this.openSource(line);
            },
        }), atom.config.onDidChange('markdown-preview-plus.markdownItConfig', () => {
            if (util_1.atomConfig().renderer === 'markdown-it')
                this.changeHandler();
        }), atom.config.onDidChange('markdown-preview-plus.pandocConfig', () => {
            if (util_1.atomConfig().renderer === 'pandoc')
                this.changeHandler();
        }), atom.config.onDidChange('markdown-preview-plus.mathConfig.latexRenderer', () => {
            util_1.handlePromise(this.handler.reload());
        }), atom.config.onDidChange('markdown-preview-plus.mathConfig.numberEquations', () => {
            util_1.handlePromise(this.handler.reload());
        }), atom.config.onDidChange('markdown-preview-plus.renderer', this.changeHandler), atom.config.onDidChange('markdown-preview-plus.useGitHubStyle', () => {
            util_1.handlePromise(this.handler.updateStyles());
        }), atom.config.onDidChange('markdown-preview-plus.syntaxThemeName', () => {
            util_1.handlePromise(this.handler.updateStyles());
        }), atom.config.onDidChange('markdown-preview-plus.importPackageStyles', () => {
            util_1.handlePromise(this.handler.updateStyles());
        }), this.handler.emitter.on('did-scroll-preview', ({ min, max }) => {
            this.didScrollPreview(min, max);
        }));
    }
    async renderMarkdown() {
        return this.renderMarkdownText(await this.getMarkdownSource());
    }
    async getHTMLToSave(savePath) {
        const source = await this.getMarkdownSource();
        return renderer.render({
            text: source,
            filePath: this.getPath(),
            grammar: this.getGrammar(),
            renderLaTeX: this.renderLaTeX,
            mode: 'save',
            savePath,
        });
    }
    async renderMarkdownText(text) {
        try {
            const domDocument = await renderer.render({
                text,
                filePath: this.getPath(),
                grammar: this.getGrammar(),
                renderLaTeX: this.renderLaTeX,
                mode: 'normal',
                imageWatcher: this.imageWatcher,
            });
            if (this.destroyed)
                return;
            await this.handler.update(domDocument.documentElement.outerHTML, this.renderLaTeX);
            await this.handler.setSourceMap(util.buildLineMap(markdownIt.getTokens(text, this.renderLaTeX)));
            this.emitter.emit('did-change-markdown');
        }
        catch (error) {
            await this.showError(error);
        }
    }
    async showError(error) {
        if (this.destroyed) {
            atom.notifications.addFatalError('Error reported on a destroyed Markdown Preview Plus view', {
                dismissable: true,
                stack: error.stack,
                detail: error.message,
            });
            return;
        }
        else if (this.loading) {
            atom.notifications.addFatalError('Error reported when Markdown Preview Plus view is loading', {
                dismissable: true,
                stack: error.stack,
                detail: error.message,
            });
            return;
        }
        else {
            return this.handler.error(error.message);
        }
    }
    async copyToClipboard() {
        await this.initialRenderPromise;
        const selection = await this.handler.getSelection();
        if (selection !== undefined) {
            atom.clipboard.write(selection);
        }
        else {
            const src = await this.getMarkdownSource();
            await util_1.copyHtml(src, this.getPath(), this.renderLaTeX);
        }
    }
}
exports.MarkdownPreviewView = MarkdownPreviewView;
MarkdownPreviewView.elementMap = new WeakMap();
//# sourceMappingURL=data:application/json;base64,eyJ2ZXJzaW9uIjozLCJmaWxlIjoibWFya2Rvd24tcHJldmlldy12aWV3LmpzIiwic291cmNlUm9vdCI6IiIsInNvdXJjZXMiOlsiLi4vLi4vc3JjL21hcmtkb3duLXByZXZpZXctdmlldy9tYXJrZG93bi1wcmV2aWV3LXZpZXcudHMiXSwibmFtZXMiOltdLCJtYXBwaW5ncyI6Ijs7QUFBQSw2QkFBNkI7QUFDN0IsK0JBQXdFO0FBQ3hFLG1DQUFpQztBQUNqQyx5QkFBeUI7QUFFekIsd0NBQXdDO0FBQ3hDLG9EQUFvRDtBQUNwRCxrQ0FBNkQ7QUFDN0QsK0JBQThCO0FBQzlCLHVEQUFrRDtBQUNsRCw4REFBb0Q7QUFDcEQsdURBQTZDO0FBQzdDLGdEQUErQztBQVEvQyxNQUFzQixtQkFBbUI7SUFrQnZDLFlBQ1UsY0FBdUIsaUJBQVUsRUFBRSxDQUFDLFVBQVU7U0FDbkQsNkJBQTZCO1FBRHhCLGdCQUFXLEdBQVgsV0FBVyxDQUNhO1FBWHhCLFlBQU8sR0FHWixJQUFJLGNBQU8sRUFBRSxDQUFBO1FBQ1IsZ0JBQVcsR0FBRyxJQUFJLDBCQUFtQixFQUFFLENBQUE7UUFDdkMsY0FBUyxHQUFHLEtBQUssQ0FBQTtRQUNuQixZQUFPLEdBQVksSUFBSSxDQUFBO1FBOEhyQixrQkFBYSxHQUFHLEdBQUcsRUFBRTtZQUM3QixvQkFBYSxDQUFDLElBQUksQ0FBQyxjQUFjLEVBQUUsQ0FBQyxDQUFBO1lBRXBDLE1BQU0sSUFBSSxHQUFHLElBQUksQ0FBQyxTQUFTLENBQUMsV0FBVyxDQUFDLElBQUksQ0FBQyxDQUFBO1lBQzdDLElBQUksSUFBSSxLQUFLLFNBQVMsSUFBSSxJQUFJLEtBQUssSUFBSSxDQUFDLFNBQVMsQ0FBQyxhQUFhLEVBQUUsRUFBRTtnQkFDakUsSUFBSSxDQUFDLFlBQVksQ0FBQyxJQUFJLENBQUMsQ0FBQTthQUN4QjtRQUNILENBQUMsQ0FBQTtRQTlIQyxJQUFJLENBQUMsb0JBQW9CLEdBQUcsSUFBSSxPQUFPLENBQUMsQ0FBQyxPQUFPLEVBQUUsRUFBRTtZQUNsRCxJQUFJLENBQUMsT0FB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ZSgncGF0aCcpXG5pbXBvcnQgeyBFbWl0dGVyLCBEaXNwb3NhYmxlLCBDb21wb3NpdGVEaXNwb3NhYmxlLCBHcmFtbWFyIH0gZnJvbSAnYXRvbSdcbmltcG9ydCB7IGRlYm91bmNlIH0gZnJvbSAnbG9kYXNoJ1xuaW1wb3J0IGZzID0gcmVxdWlyZSgnZnMnKVxuXG5pbXBvcnQgcmVuZGVyZXIgPSByZXF1aXJlKCcuLi9yZW5kZXJlcicpXG5pbXBvcnQgbWFya2Rvd25JdCA9IHJlcXVpcmUoJy4uL21hcmtkb3duLWl0LWhlbHBlcicpXG5pbXBvcnQgeyBoYW5kbGVQcm9taXNlLCBjb3B5SHRtbCwgYXRvbUNvbmZpZyB9IGZyb20gJy4uL3V0aWwnXG5pbXBvcnQgKiBhcyB1dGlsIGZyb20gJy4vdXRpbCdcbmltcG9ydCB7IFdlYnZpZXdIYW5kbGVyIH0gZnJvbSAnLi93ZWJ2aWV3LWhhbmRsZXInXG5pbXBvcnQgeyBJbWFnZVdhdGNoZXIgfSBmcm9tICcuLi9pbWFnZS13YXRjaC1oZWxwZXInXG5pbXBvcnQgeyBzYXZlQXNQREYgfSBmcm9tICcuL3BkZi1leHBvcnQtdXRpbCdcbmltcG9ydCB7IGxvYWRVc2VyTWFjcm9zIH0gZnJvbSAnLi4vbWFjcm9zLXV0aWwnXG5cbmV4cG9ydCBpbnRlcmZhY2UgU2VyaWFsaXplZE1QViB7XG4gIGRlc2VyaWFsaXplcjogJ21hcmtkb3duLXByZXZpZXctcGx1cy9NYXJrZG93blByZXZpZXdWaWV3J1xuICBlZGl0b3JJZD86IG51bWJlclxuICBmaWxlUGF0aD86IHN0cmluZ1xufVxuXG5leHBvcnQgYWJzdHJhY3QgY2xhc3MgTWFya2Rvd25QcmV2aWV3VmlldyB7XG4gIHByaXZhdGUgc3RhdGljIGVsZW1lbnRNYXAgPSBuZXcgV2Vha01hcDxIVE1MRWxlbWVudCwgTWFya2Rvd25QcmV2aWV3Vmlldz4oKVxuXG4gIHB1YmxpYyByZWFkb25seSBpbml0aWFsUmVuZGVyUHJvbWlzZTogUHJvbWlzZTx2b2lkPlxuICBwdWJsaWMgcnVuSlMhOiBNYXJrZG93blByZXZpZXdWaWV3WydoYW5kbGVyJ11bJ3J1bkpTJ11cbiAgcHJvdGVjdGVkIGhhbmRsZXIhOiBXZWJ2aWV3SGFuZGxlclxuICBwdWJsaWMgZ2V0IGVsZW1lbnQoKTogSFRNTEVsZW1lbnQge1xuICAgIHJldHVybiB0aGlzLmhhbmRsZXIuZWxlbWVudFxuICB9XG4gIHByb3RlY3RlZCBlbWl0dGVyOiBFbWl0dGVyPHtcbiAgICAnZGlkLWNoYW5nZS10aXRsZSc6IHVuZGVmaW5lZFxuICAgICdkaWQtY2hhbmdlLW1hcmtkb3duJzogdW5kZWZpbmVkXG4gIH0+ID0gbmV3IEVtaXR0ZXIoKVxuICBwcm90ZWN0ZWQgZGlzcG9zYWJsZXMgPSBuZXcgQ29tcG9zaXRlRGlzcG9zYWJsZSgpXG4gIHByb3RlY3RlZCBkZXN0cm95ZWQgPSBmYWxzZVxuICBwcml2YXRlIGxvYWRpbmc6IGJvb2xlYW4gPSB0cnVlXG4gIHByaXZhdGUgaW1hZ2VXYXRjaGVyITogSW1hZ2VXYXRjaGVyXG5cbiAgcHJvdGVjdGVkIGNvbnN0cnVjdG9yKFxuICAgIHByaXZhdGUgcmVuZGVyTGFUZVg6IGJvb2xlYW4gPSBhdG9tQ29uZmlnKCkubWF0aENvbmZpZ1xuICAgICAgLmVuYWJsZUxhdGV4UmVuZGVyaW5nQnlEZWZhdWx0LFxuICApIHtcbiAgICB0aGlzLmluaXRpYWxSZW5kZXJQcm9taXNlID0gbmV3IFByb21pc2UoKHJlc29sdmUpID0+IHtcbiAgICAgIHRoaXMuaGFuZGxlciA9IG5ldyBXZWJ2aWV3SGFuZGxlcihhc3luYyAoKSA9PiB7XG4gICAgICAgIGNvbnN0IGNvbmZpZyA9IGF0b21Db25maWcoKVxuICAgICAgICBhd2FpdCB0aGlzLmhhbmRsZXIuaW5pdCh7XG4gICAgICAgICAgdXNlck1hY3JvczogbG9hZFVzZXJNYWNyb3MoKSxcbiAgICAgICAgICBtYXRoSmF4Q29uZmlnOiBjb25maWcubWF0aENvbmZpZyxcbiAgICAgICAgICBjb250ZXh0OiAnbGl2ZS1wcmV2aWV3JyxcbiAgICAgICAgfSlcbiAgICAgICAgYXdhaXQgdGhpcy5oYW5kbGVyLnNldEJhc2VQYXRoKHRoaXMuZ2V0UGF0aCgpKVxuICAgICAgICB0aGlzLmVtaXR0ZXIuZW1pdCgnZGlkLWNoYW5nZS10aXRsZScpXG4gICAgICAgIHJlc29sdmUoXG4gICAgICAgICAgdGhpcy5yZW5kZXJNYXJrZG93bigpLnRoZW4oKCkgPT4ge1xuICAgICAgICAgICAgdGhpcy5sb2FkaW5nID0gZmFsc2VcbiAgICAgICAgICB9KSxcbiAgICAgICAgKVxuICAgICAgfSlcbiAgICAgIHRoaXMucnVuSlMgPSB0aGlzLmhhbmRsZXIucnVuSlMuYmluZCh0aGlzLmhhbmRsZXIpXG4gICAgICB0aGlzLmltYWdlV2F0Y2hlciA9IG5ldyBJbWFnZVdhdGNoZXIoXG4gICAgICAgIHRoaXMuaGFuZGxlci51cGRhdGVJbWFnZXMuYmluZCh0aGlzLmhhbmRsZXIpLFxuICAgICAgKVxuICAgICAgTWFya2Rvd25QcmV2aWV3Vmlldy5lbGVtZW50TWFwLnNldCh0aGlzLmVsZW1lbnQsIHRoaXMpXG4gICAgfSlcbiAgICB0aGlzLmhhbmRsZUV2ZW50cygpXG4gIH1cblxuICBwdWJsaWMgc3RhdGljIHZpZXdGb3JFbGVtZW50KGVsZW1lbnQ6IEhUTUxFbGVtZW50KSB7XG4gICAgcmV0dXJuIE1hcmtkb3duUHJldmlld1ZpZXcuZWxlbWVudE1hcC5nZXQoZWxlbWVudClcbiAgfVxuXG4gIHB1YmxpYyBhYnN0cmFjdCBzZXJpYWxpemUoKTogU2VyaWFsaXplZE1QVlxuXG4gIHB1YmxpYyBkZXN0cm95KCkge1xuICAgIGlmICh0aGlzLmRlc3Ryb3llZCkgcmV0dXJuXG4gICAgdGhpcy5kZXN0cm95ZWQgPSB0cnVlXG4gICAgdGhpcy5pbWFnZVdhdGNoZXIuZGlzcG9zZSgpXG4gICAgdGhpcy5kaXNwb3NhYmxlcy5kaXNwb3NlKClcbiAgICB0aGlzLmhhbmRsZXIuZGVzdHJveSgpXG4gICAgTWFya2Rvd25QcmV2aWV3Vmlldy5lbGVtZW50TWFwLmRlbGV0ZSh0aGlzLmVsZW1lbnQpXG4gIH1cblxuICBwdWJsaWMgb25EaWRDaGFuZ2VUaXRsZShjYWxsYmFjazogKCkgPT4gdm9pZCk6IERpc3Bvc2FibGUge1xuICAgIHJldHVybiB0aGlzLmVtaXR0ZXIub24oJ2RpZC1jaGFuZ2UtdGl0bGUnLCBjYWxsYmFjaylcbiAgfVxuXG4gIHB1YmxpYyBvbkRpZENoYW5nZU1hcmtkb3duKGNhbGxiYWNrOiAoKSA9PiB2b2lkKTogRGlzcG9zYWJsZSB7XG4gICAgcmV0dXJuIHRoaXMuZW1pdHRlci5vbignZGlkLWNoYW5nZS1tYXJrZG93bicsIGNhbGxiYWNrKVxuICB9XG5cbiAgcHVibGljIHRvZ2dsZVJlbmRlckxhdGV4KCkge1xuICAgIHRoaXMucmVuZGVyTGFUZVggPSAhdGhpcy5yZW5kZXJMYVRlWFxuICAgIHRoaXMuY2hhbmdlSGFuZGxlcigpXG4gIH1cblxuICBwdWJsaWMgYWJzdHJhY3QgZ2V0VGl0bGUoKTogc3RyaW5nXG5cbiAgcHVibGljIGdldERlZmF1bHRMb2NhdGlvbigpOiAnbGVmdCcgfCAncmlnaHQnIHwgJ2JvdHRvbScgfCAnY2VudGVyJyB7XG4gICAgcmV0dXJuIGF0b21Db25maWcoKS5wcmV2aWV3Q29uZmlnLnByZXZpZXdEb2NrXG4gIH1cblxuICBwdWJsaWMgZ2V0SWNvbk5hbWUoKSB7XG4gICAgcmV0dXJuICdtYXJrZG93bidcbiAgfVxuXG4gIHB1YmxpYyBhYnN0cmFjdCBnZXRVUkkoKTogc3RyaW5nXG5cbiAgcHVibGljIGFic3RyYWN0IGdldFBhdGgoKTogc3RyaW5nIHwgdW5kZWZpbmVkXG5cbiAgcHVibGljIGdldFNhdmVEaWFsb2dPcHRpb25zKCkge1xuICAgIGxldCBkZWZhdWx0UGF0aCA9IHRoaXMuZ2V0UGF0aCgpXG4gICAgaWYgKGRlZmF1bHRQYXRoID09PSB1bmRlZmluZWQpIHtcbiAgICAgIGNvbnN0IHByb2plY3RQYXRoID0gYXRvbS5wcm9qZWN0LmdldFBhdGhzKClbMF1cbiAgICAgIGRlZmF1bHRQYXRoID0gJ3VudGl0bGVkLm1kJ1xuICAgICAgaWYgKHByb2plY3RQYXRoKSB7XG4gICAgICAgIGRlZmF1bHRQYXRoID0gcGF0aC5qb2luKHByb2plY3RQYXRoLCBkZWZhdWx0UGF0aClcbiAgICAgIH1cbiAgICB9XG4gICAgZGVmYXVsdFBhdGggKz0gJy4nICsgYXRvbUNvbmZpZygpLnNhdmVDb25maWcuZGVmYXVsdFNhdmVGb3JtYXRcbiAgICByZXR1cm4geyBkZWZhdWx0UGF0aCB9XG4gIH1cblxuICBwdWJsaWMgc2F2ZUFzKGZpbGVQYXRoOiBzdHJpbmcgfCB1bmRlZmluZWQpIHtcbiAgICBpZiAoZmlsZVBhdGggPT09IHVuZGVmaW5lZCkgcmV0dXJuXG4gICAgaWYgKHRoaXMubG9hZGluZykgdGhyb3cgbmV3IEVycm9yKCdQcmV2aWV3IGlzIHN0aWxsIGxvYWRpbmcnKVxuXG4gICAgY29uc3QgeyBuYW1lLCBleHQgfSA9IHBhdGgucGFyc2UoZmlsZVBhdGgpXG5cbiAgICBpZiAoZXh0ID09PSAnLnBkZicpIHtcbiAgICAgIGhhbmRsZVByb21pc2UoXG4gICAgICAgIHRoaXMuZ2V0TWFya2Rvd25Tb3VyY2UoKS50aGVuKGFzeW5jIChtZFNvdXJjZSkgPT5cbiAgICAgICAgICBzYXZlQXNQREYoXG4gICAgICAgICAgICBtZFNvdXJjZSxcbiAgICAgICAgICAgIHRoaXMuZ2V0UGF0aCgpLFxuICAgICAgICAgICAgdGhpcy5nZXRHcmFtbWFyKCksXG4gICAgICAgICAgICB0aGlzLnJlbmRlckxhVGVYLFxuICAgICAgICAgICAgZmlsZVBhdGgsXG4gICAgICAgICAgKSxcbiAgICAgICAgKSxcbiAgICAgIClcbiAgICB9IGVsc2Uge1xuICAgICAgaGFuZGxlUHJvbWlzZShcbiAgICAgICAgdGhpcy5nZXRIVE1MVG9TYXZlKGZpbGVQYXRoKS50aGVuKGFzeW5jIChodG1sKSA9PiB7XG4gICAgICAgICAgY29uc3QgZnVsbEh0bWwgPSB1dGlsLm1rSHRtbChcbiAgICAgICAgICAgIG5hbWUsXG4gICAgICAgICAgICBodG1sLFxuICAgICAgICAgICAgdGhpcy5yZW5kZXJMYVRlWCxcbiAgICAgICAgICAgIGF3YWl0IHRoaXMuaGFuZGxlci5nZXRUZVhDb25maWcoKSxcbiAgICAgICAgICApXG5cbiAgICAgICAgICBmcy53cml0ZUZpbGVTeW5jKGZpbGVQYXRoLCBmdWxsSHRtbClcbiAgICAgICAgICByZXR1cm4gYXRvbS53b3Jrc3BhY2Uub3BlbihmaWxlUGF0aClcbiAgICAgICAgfSksXG4gICAgICApXG4gICAgfVxuICB9XG5cbiAgcHJvdGVjdGVkIGRpZFNjcm9sbFByZXZpZXcoX21pbjogbnVtYmVyLCBfbWF4OiBudW1iZXIpIHtcbiAgICAvKiBub29wLCBpbXBsZW1lbnRhdGlvbiBpbiBlZGl0b3IgcHJldmlldyAqL1xuICB9XG5cbiAgcHJvdGVjdGVkIGNoYW5nZUhhbmRsZXIgPSAoKSA9PiB7XG4gICAgaGFuZGxlUHJvbWlzZSh0aGlzLnJlbmRlck1hcmtkb3duKCkpXG5cbiAgICBjb25zdCBwYW5lID0gYXRvbS53b3Jrc3BhY2UucGFuZUZvckl0ZW0odGhpcylcbiAgICBpZiAocGFuZSAhPT0gdW5kZWZpbmVkICYmIHBhbmUgIT09IGF0b20ud29ya3NwYWNlLmdldEFjdGl2ZVBhbmUoKSkge1xuICAgICAgcGFuZS5hY3RpdmF0ZUl0ZW0odGhpcylcbiAgICB9XG4gIH1cblxuICBwcm90ZWN0ZWQgYWJzdHJhY3QgYXN5bmMgZ2V0TWFya2Rvd25Tb3VyY2UoKTogUHJvbWlzZTxzdHJpbmc+XG5cbiAgcHJvdGVjdGVkIGFic3RyYWN0IGdldEdyYW1tYXIoKTogR3JhbW1hciB8IHVuZGVmaW5lZFxuXG4gIHByb3RlY3RlZCBvcGVuU291cmNlKGluaXRpYWxMaW5lPzogbnVtYmVyKSB7XG4gICAgY29uc3QgcGF0aCA9IHRoaXMuZ2V0UGF0aCgpXG4gICAgaWYgKHBhdGggPT09IHVuZGVmaW5lZCkgcmV0dXJuXG4gICAgaGFuZGxlUHJvbWlzZShcbiAgICAgIGF0b20ud29ya3NwYWNlLm9wZW4ocGF0aCwge1xuICAgICAgICBpbml0aWFsTGluZSxcbiAgICAgICAgc2VhcmNoQWxsUGFuZXM6IHRydWUsXG4gICAgICB9KSxcbiAgICApXG4gIH1cblxuICBwcm90ZWN0ZWQgc3luY1ByZXZpZXcobGluZTogbnVtYmVyLCBmbGFzaDogYm9vbGVhbikge1xuICAgIGhhbmRsZVByb21pc2UodGhpcy5oYW5kbGVyLnN5bmMobGluZSwgZmxhc2gpKVxuICB9XG5cbiAgcHJvdGVjdGVkIG9wZW5OZXdXaW5kb3coKSB7XG4gICAgY29uc3QgcGF0aCA9IHRoaXMuZ2V0UGF0aCgpXG4gICAgaWYgKCFwYXRoKSB7XG4gICAgICBhdG9tLm5vdGlmaWNhdGlvbnMuYWRkV2FybmluZyhcbiAgICAgICAgJ0NhbiBub3Qgb3BlbiB0aGlzIHByZXZpZXcgaW4gbmV3IHdpbmRvdzogbm8gZmlsZSBwYXRoJyxcbiAgICAgIClcbiAgICAgIHJldHVyblxuICAgIH1cbiAgICBhdG9tLm9wZW4oe1xuICAgICAgcGF0aHNUb09wZW46IFtgbWFya2Rvd24tcHJldmlldy1wbHVzOi8vZmlsZS8ke3BhdGh9YF0sXG4gICAgICBuZXdXaW5kb3c6IHRydWUsXG4gICAgfSlcbiAgICB1dGlsLmRlc3Ryb3kodGhpcylcbiAgfVxuXG4gIHByaXZhdGUgaGFuZGxlRXZlbnRzKCkge1xuICAgIHRoaXMuZGlzcG9zYWJsZXMuYWRkKFxuICAgICAgLy8gYXRvbSBldmVudHNcbiAgICAgIGF0b20uZ3JhbW1hcnMub25EaWRBZGRHcmFtbWFyKCgpID0+XG4gICAgICAgIGRlYm91bmNlKCgpID0+IHtcbiAgICAgICAgICBoYW5kbGVQcm9taXNlKHRoaXMucmVuZGVyTWFya2Rvd24oKSlcbiAgICAgICAgfSwgMjUwKSxcbiAgICAgICksXG4gICAgICBhdG9tLmdyYW1tYXJzLm9uRGlkVXBkYXRlR3JhbW1hcihcbiAgICAgICAgZGVib3VuY2UoKCkgPT4ge1xuICAgICAgICAgIGhhbmRsZVByb21pc2UodGhpcy5yZW5kZXJNYXJrZG93bigpKVxuICAgICAgICB9LCAyNTApLFxuICAgICAgKSxcbiAgICAgIGF0b20uY29tbWFuZHMuYWRkKHRoaXMuZWxlbWVudCwge1xuICAgICAgICAnY29yZTptb3ZlLXVwJzogKCkgPT4gdGhpcy5lbGVtZW50LnNjcm9sbEJ5KHsgdG9wOiAtMTAgfSksXG4gICAgICAgICdjb3JlOm1vdmUtZG93bic6ICgpID0+IHRoaXMuZWxlbWVudC5zY3JvbGxCeSh7IHRvcDogMTAgfSksXG4gICAgICAgICdjb3JlOmNvcHknOiAoKSA9PiB7XG4gICAgICAgICAgaGFuZGxlUHJvbWlzZSh0aGlzLmNvcHlUb0NsaXBib2FyZCgpKVxuICAgICAgICB9LFxuICAgICAgICAnbWFya2Rvd24tcHJldmlldy1wbHVzOm9wZW4tZGV2LXRvb2xzJzogKCkgPT4ge1xuICAgICAgICAgIHRoaXMuaGFuZGxlci5vcGVuRGV2VG9vbHMoKVxuICAgICAgICB9LFxuICAgICAgICAnbWFya2Rvd24tcHJldmlldy1wbHVzOm5ldy13aW5kb3cnOiAoKSA9PiB7XG4gICAgICAgICAgdGhpcy5vcGVuTmV3V2luZG93KClcbiAgICAgICAgfSxcbiAgICAgICAgJ21hcmtkb3duLXByZXZpZXctcGx1czpwcmludCc6ICgpID0+IHtcbiAgICAgICAgICB0aGlzLmhhbmRsZXIucHJpbnQoKVxuICAgICAgICB9LFxuICAgICAgICAnbWFya2Rvd24tcHJldmlldy1wbHVzOnpvb20taW4nOiAoKSA9PiB7XG4gICAgICAgICAgdGhpcy5oYW5kbGVyLnpvb21JbigpXG4gICAgICAgIH0sXG4gICAgICAgICdtYXJrZG93bi1wcmV2aWV3LXBsdXM6em9vbS1vdXQnOiAoKSA9PiB7XG4gICAgICAgICAgdGhpcy5oYW5kbGVyLnpvb21PdXQoKVxuICAgICAgICB9LFxuICAgICAgICAnbWFya2Rvd24tcHJldmlldy1wbHVzOnJlc2V0LXpvb20nOiAoKSA9PiB7XG4gICAgICAgICAgdGhpcy5oYW5kbGVyLnJlc2V0Wm9vbSgpXG4gICAgICAgIH0sXG4gICAgICAgICdtYXJrZG93bi1wcmV2aWV3LXBsdXM6c3luYy1zb3VyY2UnOiBhc3luYyAoX2V2ZW50KSA9PiB7XG4gICAgICAgICAgY29uc3QgbGluZSA9IGF3YWl0IHRoaXMuaGFuZGxlci5zeW5jU291cmNlKClcbiAgICAgICAgICB0aGlzLm9wZW5Tb3VyY2UobGluZSlcbiAgICAgICAgfSxcbiAgICAgIH0pLFxuICAgICAgYXRvbS5jb25maWcub25EaWRDaGFuZ2UoJ21hcmtkb3duLXByZXZpZXctcGx1cy5tYXJrZG93bkl0Q29uZmlnJywgKCkgPT4ge1xuICAgICAgICBpZiAoYXRvbUNvbmZpZygpLnJlbmRlcmVyID09PSAnbWFya2Rvd24taXQnKSB0aGlzLmNoYW5nZUhhbmRsZXIoKVxuICAgICAgfSksXG4gICAgICBhdG9tLmNvbmZpZy5vbkRpZENoYW5nZSgnbWFya2Rvd24tcHJldmlldy1wbHVzLnBhbmRvY0NvbmZpZycsICgpID0+IHtcbiAgICAgICAgaWYgKGF0b21Db25maWcoKS5yZW5kZXJlciA9PT0gJ3BhbmRvYycpIHRoaXMuY2hhbmdlSGFuZGxlcigpXG4gICAgICB9KSxcbiAgICAgIGF0b20uY29uZmlnLm9uRGlkQ2hhbmdlKFxuICAgICAgICAnbWFya2Rvd24tcHJldmlldy1wbHVzLm1hdGhDb25maWcubGF0ZXhSZW5kZXJlcicsXG4gICAgICAgICgpID0+IHtcbiAgICAgICAgICBoYW5kbGVQcm9taXNlKHRoaXMuaGFuZGxlci5yZWxvYWQoKSlcbiAgICAgICAgfSxcbiAgICAgICksXG4gICAgICBhdG9tLmNvbmZpZy5vbkRpZENoYW5nZShcbiAgICAgICAgJ21hcmtkb3duLXByZXZpZXctcGx1cy5tYXRoQ29uZmlnLm51bWJlckVxdWF0aW9ucycsXG4gICAgICAgICgpID0+IHtcbiAgICAgICAgICBoYW5kbGVQcm9taXNlKHRoaXMuaGFuZGxlci5yZWxvYWQoKSlcbiAgICAgICAgfSxcbiAgICAgICksXG4gICAgICBhdG9tLmNvbmZpZy5vbkRpZENoYW5nZShcbiAgICAgICAgJ21hcmtkb3duLXByZXZpZXctcGx1cy5yZW5kZXJlcicsXG4gICAgICAgIHRoaXMuY2hhbmdlSGFuZGxlcixcbiAgICAgICksXG4gICAgICBhdG9tLmNvbmZpZy5vbkRpZENoYW5nZSgnbWFya2Rvd24tcHJldmlldy1wbHVzLnVzZUdpdEh1YlN0eWxlJywgKCkgPT4ge1xuICAgICAgICBoYW5kbGVQcm9taXNlKHRoaXMuaGFuZGxlci51cGRhdGVTdHlsZXMoKSlcbiAgICAgIH0pLFxuICAgICAgYXRvbS5jb25maWcub25EaWRDaGFuZ2UoJ21hcmtkb3duLXByZXZpZXctcGx1cy5zeW50YXhUaGVtZU5hbWUnLCAoKSA9PiB7XG4gICAgICAgIGhhbmRsZVByb21pc2UodGhpcy5oYW5kbGVyLnVwZGF0ZVN0eWxlcygpKVxuICAgICAgfSksXG4gICAgICBhdG9tLmNvbmZpZy5vbkRpZENoYW5nZShcbiAgICAgICAgJ21hcmtkb3duLXByZXZpZXctcGx1cy5pbXBvcnRQYWNrYWdlU3R5bGVzJyxcbiAgICAgICAgKCkgPT4ge1xuICAgICAgICAgIGhhbmRsZVByb21pc2UodGhpcy5oYW5kbGVyLnVwZGF0ZVN0eWxlcygpKVxuICAgICAgICB9LFxuICAgICAgKSxcblxuICAgICAgLy8gd2VidmlldyBldmVudHNcbiAgICAgIHRoaXMuaGFuZGxlci5lbWl0dGVyLm9uKCdkaWQtc2Nyb2xsLXByZXZpZXcnLCAoeyBtaW4sIG1heCB9KSA9PiB7XG4gICAgICAgIHRoaXMuZGlkU2Nyb2xsUHJldmlldyhtaW4sIG1heClcbiAgICAgIH0pLFxuICAgIClcbiAgfVxuXG4gIHByaXZhdGUgYXN5bmMgcmVuZGVyTWFya2Rvd24oKTogUHJvbWlzZTx2b2lkPiB7XG4gICAgcmV0dXJuIHRoaXMucmVuZGVyTWFya2Rvd25UZXh0KGF3YWl0IHRoaXMuZ2V0TWFya2Rvd25Tb3VyY2UoKSlcbiAgfVxuXG4gIHByaXZhdGUgYXN5bmMgZ2V0SFRNTFRvU2F2ZShzYXZlUGF0aDogc3RyaW5nKSB7XG4gICAgY29uc3Qgc291cmNlID0gYXdhaXQgdGhpcy5nZXRNYXJrZG93blNvdXJjZSgpXG4gICAgcmV0dXJuIHJlbmRlcmVyLnJlbmRlcih7XG4gICAgICB0ZXh0OiBzb3VyY2UsXG4gICAgICBmaWxlUGF0aDogdGhpcy5nZXRQYXRoKCksXG4gICAgICBncmFtbWFyOiB0aGlzLmdldEdyYW1tYXIoKSxcbiAgICAgIHJlbmRlckxhVGVYOiB0aGlzLnJlbmRlckxhVGVYLFxuICAgICAgbW9kZTogJ3NhdmUnLFxuICAgICAgc2F2ZVBhdGgsXG4gICAgfSlcbiAgfVxuXG4gIHByaXZhdGUgYXN5bmMgcmVuZGVyTWFya2Rvd25UZXh0KHRleHQ6IHN0cmluZyk6IFByb21pc2U8dm9pZD4ge1xuICAgIHRyeSB7XG4gICAgICBjb25zdCBkb21Eb2N1bWVudCA9IGF3YWl0IHJlbmRlcmVyLnJlbmRlcih7XG4gICAgICAgIHRleHQsXG4gICAgICAgIGZpbGVQYXRoOiB0aGlzLmdldFBhdGgoKSxcbiAgICAgICAgZ3JhbW1hcjogdGhpcy5nZXRHcmFtbWFyKCksXG4gICAgICAgIHJlbmRlckxhVGVYOiB0aGlzLnJlbmRlckxhVGVYLFxuICAgICAgICBtb2RlOiAnbm9ybWFsJyxcbiAgICAgICAgaW1hZ2VXYXRjaGVyOiB0aGlzLmltYWdlV2F0Y2hlcixcbiAgICAgIH0pXG5cbiAgICAgIGlmICh0aGlzLmRlc3Ryb3llZCkgcmV0dXJuXG4gICAgICBhd2FpdCB0aGlzLmhhbmRsZXIudXBkYXRlKFxuICAgICAgICBkb21Eb2N1bWVudC5kb2N1bWVudEVsZW1lbnQhLm91dGVySFRNTCxcbiAgICAgICAgdGhpcy5yZW5kZXJMYVRlWCxcbiAgICAgIClcbiAgICAgIGF3YWl0IHRoaXMuaGFuZGxlci5zZXRTb3VyY2VNYXAoXG4gICAgICAgIHV0aWwuYnVpbGRMaW5lTWFwKG1hcmtkb3duSXQuZ2V0VG9rZW5zKHRleHQsIHRoaXMucmVuZGVyTGFUZVgpKSxcbiAgICAgIClcbiAgICAgIHRoaXMuZW1pdHRlci5lbWl0KCdkaWQtY2hhbmdlLW1hcmtkb3duJylcbiAgICB9IGNhdGNoIChlcnJvcikge1xuICAgICAgYXdhaXQgdGhpcy5zaG93RXJyb3IoZXJyb3IgYXMgRXJyb3IpXG4gICAgfVxuICB9XG5cbiAgcHJpdmF0ZSBhc3luYyBzaG93RXJyb3IoZXJyb3I6IEVycm9yKSB7XG4gICAgaWYgKHRoaXMuZGVzdHJveWVkKSB7XG4gICAgICBhdG9tLm5vdGlmaWNhdGlvbnMuYWRkRmF0YWxFcnJvcihcbiAgICAgICAgJ0Vycm9yIHJlcG9ydGVkIG9uIGEgZGVzdHJveWVkIE1hcmtkb3duIFByZXZpZXcgUGx1cyB2aWV3JyxcbiAgICAgICAge1xuICAgICAgICAgIGRpc21pc3NhYmxlOiB0cnVlLFxuICAgICAgICAgIHN0YWNrOiBlcnJvci5zdGFjayxcbiAgICAgICAgICBkZXRhaWw6IGVycm9yLm1lc3NhZ2UsXG4gICAgICAgIH0sXG4gICAgICApXG4gICAgICByZXR1cm5cbiAgICB9IGVsc2UgaWYgKHRoaXMubG9hZGluZykge1xuICAgICAgYXRvbS5ub3RpZmljYXRpb25zLmFkZEZhdGFsRXJyb3IoXG4gICAgICAgICdFcnJvciByZXBvcnRlZCB3aGVuIE1hcmtkb3duIFByZXZpZXcgUGx1cyB2aWV3IGlzIGxvYWRpbmcnLFxuICAgICAgICB7XG4gICAgICAgICAgZGlzbWlzc2FibGU6IHRydWUsXG4gICAgICAgICAgc3RhY2s6IGVycm9yLnN0YWNrLFxuICAgICAgICAgIGRldGFpbDogZXJyb3IubWVzc2FnZSxcbiAgICAgICAgfSxcbiAgICAgIClcbiAgICAgIHJldHVyblxuICAgIH0gZWxzZSB7XG4gICAgICByZXR1cm4gdGhpcy5oYW5kbGVyLmVycm9yKGVycm9yLm1lc3NhZ2UpXG4gICAgfVxuICB9XG5cbiAgcHJpdmF0ZSBhc3luYyBjb3B5VG9DbGlwYm9hcmQoKTogUHJvbWlzZTx2b2lkPiB7XG4gICAgYXdhaXQgdGhpcy5pbml0aWFsUmVuZGVyUHJvbWlzZVxuICAgIGNvbnN0IHNlbGVjdGlvbiA9IGF3YWl0IHRoaXMuaGFuZGxlci5nZXRTZWxlY3Rpb24oKVxuICAgIC8vIFVzZSBzdHVwaWQgY29weSBldmVudCBoYW5kbGVyIGlmIHRoZXJlIGlzIHNlbGVjdGVkIHRleHQgaW5zaWRlIHRoaXMgdmlld1xuICAgIGlmIChzZWxlY3Rpb24gIT09IHVuZGVmaW5lZCkge1xuICAgICAgLy8gVE9ETzogcmljaCBjbGlwYm9hcmQgc3VwcG9ydFxuICAgICAgYXRvbS5jbGlwYm9hcmQud3JpdGUoc2VsZWN0aW9uKVxuICAgIH0gZWxzZSB7XG4gICAgICBjb25zdCBzcmMgPSBhd2FpdCB0aGlzLmdldE1hcmtkb3duU291cmNlKClcbiAgICAgIGF3YWl0IGNvcHlIdG1sKHNyYywgdGhpcy5nZXRQYXRoKCksIHRoaXMucmVuZGVyTGFUZVgpXG4gICAgfVxuICB9XG59XG4iXX0=